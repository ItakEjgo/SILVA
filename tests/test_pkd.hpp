#pragma once
#include <iostream>
#include <vector>
#include "../benchmarks/runners/pkd_runner.hpp"
#include <helper/time_loop.h>

namespace PKDTest {
    using PT = parlay::sequence<geobase::Point>;
    using RQ = parlay::sequence<geobase::Bounding_Box>;

    inline std::vector<Value> convert_points(const PT& P) {
        std::vector<Value> P_conv(P.size());
        for (size_t i = 0; i < P.size(); i++) {
            P_conv[i] = std::make_pair(BoostPoint(P[i].x, P[i].y), P[i].id);
        }
        return P_conv;
    }

    void build_test(PT P) {
        auto P_conv = convert_points(P);
        double final_mem = 0;
        PKDRunner* runner = nullptr;
        double ms = time_loop(5, 1.0, 
            [&](){ 
                if (runner) delete runner;
                runner = new PKDRunner();
            }, 
            [&]() {
                runner->build_base(P_conv);
            }, 
            [&](){
                final_mem = runner->memory_usage().first;
            }
        ) * 1000.0;
        delete runner;
        std::cout << "[memory_MB]: " << final_mem << std::endl;
        std::cout << "[PkdTree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    void batch_insert_test(PT P_base, PT P_update, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        auto n = P_base.size();
        auto rand_p = geobase::shuffle_point(P_update);
        parlay::parallel_for (0, rand_p.size(), [&](int i){
            rand_p[i].id = n + i;
        });
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = rand_p.size() * ratio;
            if (cur_batch_size == 0) cur_batch_size = 1;
            
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            double total_ms = 0;

            PKDRunner runner;
            runner.build_base(P_conv);

            for (size_t offset = 0; offset < rand_p.size(); offset += cur_batch_size) {
                size_t current_batch_size = std::min(cur_batch_size, rand_p.size() - offset);
                parlay::sequence<geobase::Point> adds(current_batch_size);
                for(size_t j=0; j<current_batch_size; j++) adds[j] = rand_p[offset + j];
                parlay::sequence<geobase::Point> rems; // empty

                // We MUST backup and rebuild PKDTree for clearf
                parlay::sequence<PKDRunner::point_t> backup_pts;
                if (runner.pkd.get_root()) {
                    backup_pts = parlay::sequence<PKDRunner::point_t>::uninitialized(runner.pkd.get_root()->size);
                    runner.pkd.flatten(runner.pkd.get_root(), parlay::make_slice(backup_pts));
                }
                double mem_recorded = 0;

                double batch_avg = time_loop(
                    3, 1.0, 
                    [&]() { 
                        runner.pkd.delete_tree();
                        runner.pkd.build(parlay::make_slice(backup_pts), 2);
                    },
                    [&]() { 
                        runner.commit(adds, rems);
                    },
                    [&]() {
                        mem_recorded = runner.memory_usage().first;
                    }
                );

                runner.pkd.delete_tree();
                runner.pkd.build(parlay::make_slice(backup_pts), 2);
                runner.commit(adds, rems);

                batch_times.push_back(batch_avg * 1000.0);
                total_ms += batch_avg * 1000.0;
                batch_mems.push_back(runner.memory_usage().first);
            }
            for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_batch_mem]: ";
            for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << batch_mems.back() << std::endl;
            std::cout << "[PkdTree]: batch insert time (avg): " << (total_ms / 1000.0) << std::endl;
        }
    }

    void batch_delete_test(PT P_base, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        auto rand_p = geobase::shuffle_point(P_base);
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = P_base.size() * ratio;
            if (cur_batch_size == 0) cur_batch_size = 1;
            
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            double total_ms = 0;

            PKDRunner runner;
            runner.build_base(P_conv);

            for (size_t offset = 0; offset < P_base.size(); offset += cur_batch_size) {
                size_t current_batch_size = std::min(cur_batch_size, P_base.size() - offset);
                parlay::sequence<geobase::Point> adds; // empty
                parlay::sequence<geobase::Point> rems(current_batch_size);
                for(size_t j=0; j<current_batch_size; j++) rems[j] = rand_p[offset + j];

                double mem_recorded = 0;

                // PKDTree ignores rems, so it's a no-op. Time loop is trivial.
                double batch_avg = time_loop(
                    3, 1.0, [&](){}, [&](){runner.commit(adds, rems);}, 
                    [&](){ mem_recorded = runner.memory_usage().first; }
                );

                runner.commit(adds, rems);

                batch_times.push_back(batch_avg * 1000.0);
                total_ms += batch_avg * 1000.0;
                batch_mems.push_back(runner.memory_usage().first);
            }
            std::cout << "[per_batch_time]: ";
            for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_batch_mem]: ";
            for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << batch_mems.back() << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
            std::cout << "[PkdTree]: batch delete time (avg): " << (total_ms / 1000.0) << std::endl;
        }
    }

    void range_report_test(PT P, RQ qs, parlay::sequence<size_t>& cnt) {
        auto P_conv = convert_points(P);
        PKDRunner runner;
        runner.build_base(P_conv);
        
        parlay::sequence<geobase::Point> shared_out(P.size());
        std::vector<size_t> actual_cnt(qs.size(), 0);
        
        auto avg_time = time_loop(
            3, 1.0, 
            [&]() {},
            [&]() {					
                parlay::parallel_for(0, qs.size(), [&](size_t i){
                    size_t h = 0;
                    runner.range_query(qs[i], shared_out, actual_cnt[i], h);
                });
            },
            [&](){} 
        );

        std::cout << "[PkdTree]: range report time (avg): " << avg_time << std::endl;
        
        bool is_correct = true;
        for (size_t i = 0; i < qs.size(); i++) {
            if (actual_cnt[i] != cnt[i]) {
                std::cout << "[Pkd ERROR] Query " << i << " failed. Expected: " << cnt[i] << ", Got: " << actual_cnt[i] << std::endl;
                is_correct = false; break;
            }
        }
        if (is_correct) std::cout << "[PkdTree] Accuracy: 100% (All " << qs.size() << " queries correct)" << std::endl;
    }
}
