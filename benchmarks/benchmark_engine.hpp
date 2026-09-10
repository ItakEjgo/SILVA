#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <parlay/sequence.h>
#include <silva/geo/point.hpp>
#include <silva/geo/operations.hpp>
#include <silva/core/dataset_loader.hpp>
#include <sys/stat.h>
#include <cpam/get_time.h>

// The BenchmarkEngine takes a DatasetLoader (which contains all events in memory)
// and replays the events, maintaining active_nodes and all_ever_existed_nodes.
// Every `q_step` commits, it pauses, generates queries, and calls the Runner to test them.

template <typename Runner>
class BenchmarkEngine {
public:
    static void run_workload(Runner& runner, DatasetLoader& data, const std::string& run_algo, size_t q_step) {
        std::unordered_map<size_t, geobase::Point> active_nodes;
        std::vector<geobase::Point> all_ever_existed_nodes;
        
        // 1. Initialize Base Data
        for(auto& pt : data.base_data) {
            active_nodes[pt.id] = pt;
            all_ever_existed_nodes.push_back(pt);
        }
        
        std::cout << "[Event-Driven] Executing Workload for " << Runner::name() << " with step size: " << q_step << std::endl;
        
        // Ensure reproducible random queries for EVERY runner
        srand(42); 
        
        size_t commit_count = 0;
        double accumulated_commit_ms = 0.0;
        
        // Ensure result directory exists
        mkdir("verification_results", 0777);
        std::ofstream commit_fout("verification_results/CommitLog_" + run_algo + ".txt");
        commit_fout << "Algo | CS_ID | Year | Adds | Rems | Time_ms | Mem_Tree_MB | Mem_Delta_MB | DeltaData_MB\n";

        for (const auto& event : data.all_events) {
            std::vector<geobase::Point> adds;
            std::vector<geobase::Point> rems;
            
            for (const auto& op : event.ops) {
                if (op.op_type == 'A') {
                    geobase::Point pt(op.node_id, op.lon, op.lat);
                    adds.push_back(pt);
                    active_nodes[pt.id] = pt;
                    all_ever_existed_nodes.push_back(pt);
                } else if (op.op_type == 'M') {
                    if (active_nodes.find(op.node_id) != active_nodes.end()) {
                        rems.push_back(active_nodes[op.node_id]);
                    }
                    geobase::Point pt(op.node_id, op.lon, op.lat);
                    adds.push_back(pt);
                    active_nodes[pt.id] = pt;
                    all_ever_existed_nodes.push_back(pt);
                } else if (op.op_type == 'D') {
                    if (active_nodes.find(op.node_id) != active_nodes.end()) {
                        rems.push_back(active_nodes[op.node_id]);
                        active_nodes.erase(op.node_id);
                    }
                }
            }
            
            parlay::sequence<geobase::Point> pre_sort_adds(adds.begin(), adds.end());
            parlay::sequence<geobase::Point> pre_sort_rems(rems.begin(), rems.end());
            auto sorted_adds_vec = geobase::get_sorted_points(pre_sort_adds);
            auto sorted_rems_vec = geobase::get_sorted_points(pre_sort_rems);
            
            parlay::sequence<geobase::Point> seq_adds(sorted_adds_vec.begin(), sorted_adds_vec.end());
            parlay::sequence<geobase::Point> seq_rems(sorted_rems_vec.begin(), sorted_rems_vec.end());

            // Commit the runner
            cpam::timer t_c;
            runner.commit(seq_adds, seq_rems);
            double cur_ms = t_c.stop() * 1000.0;
            accumulated_commit_ms += cur_ms;
            
            double delta_data_mb = (adds.size() + rems.size()) * sizeof(geobase::Point) / (1024.0 * 1024.0);
            auto mems = runner.memory_usage();
            commit_fout << Runner::name() << " | " << event.changeset_id << " | " << event.year << " | " 
                        << adds.size() << " | " << rems.size() << " | " << cur_ms << " | " 
                        << mems.first << " | 0.0 | " << delta_data_mb << "\n";
            
            commit_count++;
            
            // Execute checkpoint queries
            if (commit_count > 0 && commit_count % q_step == 0) {
                std::cout << "[Event-Driven] Executing Checkpoint " << commit_count << "..." << std::endl;
                std::string fname = "verification_results/Checkpoint_" + std::to_string(commit_count) + "_" + run_algo + ".txt";
                
                std::ofstream fout(fname);
                fout << "Algo | QType | QID | Count | Hash | Time_ms | Mem_Tree_MB | Mem_Delta_MB | Nodes\n";
                fout.close();
                
                run_hybrid_queries(runner, fname, active_nodes, all_ever_existed_nodes, data.largest_mbr);
                
                accumulated_commit_ms = 0.0;
            }
        }
        
        std::ofstream end_flag("verification_results/2026_EndQuery_" + run_algo + ".txt");
        end_flag << "Completed Workload." << std::endl;
    }

private:
    static void run_hybrid_queries(Runner& runner, const std::string& fname, 
                                   std::unordered_map<size_t, geobase::Point>& active_nodes, 
                                   std::vector<geobase::Point>& all_ever_existed_nodes,
                                   const geobase::Bounding_Box& world_mbr) {
        std::ofstream fout(fname, std::ios_base::app);
        if (!fout.is_open()) return;
        
        size_t n = active_nodes.size();
        if (n == 0) return;
        
        // Dynamically pre-allocate result buffer to prevent OOM
        parlay::sequence<geobase::Point> shared_out(n + 100000);
        
        parlay::sequence<geobase::Point> current_knn(100);
        for(int i=0; i<100; i++) {
            while(true) {
                auto pt = all_ever_existed_nodes[rand() % all_ever_existed_nodes.size()];
                if (active_nodes.find(pt.id) != active_nodes.end()) {
                    current_knn[i] = active_nodes[pt.id];
                    break;
                }
            }
        }
        
        parlay::sequence<geobase::Bounding_Box> qs_small(100), qs_med(100), qs_large(100);
        
        double world_w = world_mbr.second.x - world_mbr.first.x;
        double world_h = world_mbr.second.y - world_mbr.first.y;
        
        auto build_qs = [&](double ratio, parlay::sequence<geobase::Bounding_Box>& qs) {
            double dx = world_w * ratio;
            double dy = world_h * ratio;
            for(int i=0; i<100; i++) {
                geobase::Point center;
                while(true) {
                    auto pt = all_ever_existed_nodes[rand() % all_ever_existed_nodes.size()];
                    if (active_nodes.find(pt.id) != active_nodes.end()) {
                        center = active_nodes[pt.id];
                        break;
                    }
                }
                qs[i] = geobase::Bounding_Box(
                    geobase::Point(center.x - dx, center.y - dy),
                    geobase::Point(center.x + dx, center.y + dy)
                );
            }
        };
        build_qs(0.0001, qs_small);
        build_qs(0.001, qs_med);
        build_qs(0.01, qs_large);
        
        auto run_range_set = [&](const parlay::sequence<geobase::Bounding_Box>& qs, std::string q_label) {
            for(size_t i=0; i<qs.size(); i++) {
                geobase::Bounding_Box q_copy = qs[i];
                size_t cnt = 0;
                size_t h = 0;
                cpam::timer t_q;
                for (int rep = 0; rep < 3; rep++) {
                    runner.reset_query_stats();
                    runner.range_query(q_copy, shared_out, cnt, h);
                }
                double q_ms = (t_q.stop() * 1000.0) / 3.0;
                auto mems = runner.memory_usage();
                fout << Runner::name() << " | " << q_label << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << runner.get_nodes_touched() << "\n";
            }
        };
        run_range_set(qs_small, "Range_Small");
        run_range_set(qs_med, "Range_Medium");
        run_range_set(qs_large, "Range_Large");
        
        auto run_knn_set = [&](int K, std::string qk) {
            for(size_t i=0; i<100; i++) {
                size_t cnt = 0;
                size_t h = 0;
                cpam::timer t_q;
                for (int rep = 0; rep < 3; rep++) {
                    runner.reset_query_stats();
                    runner.knn_query(current_knn[i], K, cnt, h);
                }
                double q_ms = (t_q.stop() * 1000.0) / 3.0;
                auto mems = runner.memory_usage();
                fout << Runner::name() << " | " << qk << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << runner.get_nodes_touched() << "\n";
            }
        };
        run_knn_set(10, "kNN_10");
        run_knn_set(100, "kNN_100");
    }
};
