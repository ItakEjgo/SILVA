#pragma once
#include <iostream>
#include <vector>
#include <silva/baselines/pkdlog_tree.hpp>
#include <helper/time_loop.h>

namespace PKDLogTest {
    using PT = parlay::sequence<geobase::Point>;
    using RQ = parlay::sequence<geobase::Bounding_Box>;

    inline std::vector<Value> convert_points(const PT& P) {
        std::vector<Value> P_conv(P.size());
        for (size_t i = 0; i < P.size(); i++) {
            P_conv[i] = std::make_pair(BoostPoint(P[i].x, P[i].y), P[i].id);
        }
        return P_conv;
    }

    inline void build_test(PT P) {
        auto P_conv = convert_points(P);
        double final_mem = 0;
        std::shared_ptr<PKDLog::VersionNode> tree;
        double ms = time_loop(5, 1.0, 
            [&](){ tree.reset(); },
            [&]() {
                tree = PKDLog::map_init(P_conv);
            }, 
            [&](){
                final_mem = cpam::cpam_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
            }
        ) * 1000.0;
        std::cout << "[memory_MB]: " << final_mem << std::endl;
        std::cout << "[PkdLogTree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    inline void batch_insert_test(PT P_base, PT P_update, parlay::sequence<double>& batch_ratios, double p = 1.0) {
        auto P_conv = convert_points(P_base);
        auto n = P_base.size();
        auto rand_p = geobase::shuffle_point(P_update);
        parlay::parallel_for (0, rand_p.size(), [&](int i){
            rand_p[i].id = n + i;
        });
        auto P_update_conv = convert_points(rand_p);
        
        auto base_ver = PKDLog::map_init(P_conv);
        
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = std::max<size_t>(1, P_update_conv.size() * ratio);
            
            std::vector<std::shared_ptr<PKDLog::VersionNode>> versions;
            versions.push_back(base_ver);
            
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            double total_ms = 0;
            
            cout << "[Testing Ratio]: " << ratio << endl;
            for (size_t offset = 0; offset < P_update_conv.size(); offset += cur_batch_size) {
                size_t current_batch_size = std::min(cur_batch_size, P_update_conv.size() - offset);
                std::vector<Value> batch(P_update_conv.begin() + offset, P_update_conv.begin() + offset + current_batch_size);
                
                std::shared_ptr<PKDLog::VersionNode> test_ver;
                
                double batch_avg = time_loop(
                    3, 1.0, 
                    [&]() { test_ver.reset(); },
                    [&]() { test_ver = PKDLog::map_insert(versions.back(), batch, p); },
                    [&]() {}
                );
                
                batch_times.push_back(batch_avg * 1000.0);
                total_ms += batch_avg * 1000.0;
                versions.push_back(test_ver);
                double mem_mb = cpam::cpam_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
                batch_mems.push_back(mem_mb);
                cout << "[step_time]: " << batch_avg * 1000.0 << " [step_mem]: " << mem_mb << endl;
            }
            
            std::cout << "[per_batch_time]: ";
            for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_batch_mem]: ";
            for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << batch_mems.back() << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
            std::cout << "[PkdLogTree]: batch insert time (avg): " << (total_ms / 1000.0) << std::endl;
        }
    }

    inline void batch_delete_test(PT P_base, parlay::sequence<double>& batch_ratios, double p = 1.0) {
        auto rand_p = geobase::shuffle_point(P_base);
        auto P_conv = convert_points(P_base);
        auto P_delete_conv = convert_points(rand_p);
        
        auto base_ver = PKDLog::map_init(P_conv);
        
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = std::max<size_t>(1, P_base.size() * ratio);
            
            std::vector<std::shared_ptr<PKDLog::VersionNode>> versions;
            versions.push_back(base_ver);
            
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            double total_ms = 0;
            
            cout << "[Testing Ratio]: " << ratio << endl;
            for (size_t offset = 0; offset < P_delete_conv.size(); offset += cur_batch_size) {
                size_t current_batch_size = std::min(cur_batch_size, P_delete_conv.size() - offset);
                std::vector<Value> batch(P_delete_conv.begin() + offset, P_delete_conv.begin() + offset + current_batch_size);
                
                std::shared_ptr<PKDLog::VersionNode> test_ver;
                
                double batch_avg = time_loop(
                    3, 1.0, 
                    [&]() { test_ver.reset(); },
                    [&]() { test_ver = PKDLog::map_delete(versions.back(), batch, p); },
                    [&]() {}
                );
                
                batch_times.push_back(batch_avg * 1000.0);
                total_ms += batch_avg * 1000.0;
                versions.push_back(test_ver);
                double mem_mb = cpam::cpam_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
                batch_mems.push_back(mem_mb);
                cout << "[step_time]: " << batch_avg * 1000.0 << " [step_mem]: " << mem_mb << endl;
            }
            
            std::cout << "[per_batch_time]: ";
            for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_batch_mem]: ";
            for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << batch_mems.back() << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
            std::cout << "[PkdLogTree]: batch delete time (avg): " << (total_ms / 1000.0) << std::endl;
        }
    }

    inline void range_report_test(PT P, RQ qs, parlay::sequence<size_t>& cnt) {
        auto P_conv = convert_points(P);
        auto tree = PKDLog::map_init(P_conv);
        double ms = 0;
        
        std::vector<size_t> actual_cnt(qs.size(), 0);
        for (int rep = 0; rep < 3; rep++) {
            parlay::internal::timer t;
            for (size_t i = 0; i < qs.size(); i++) {
                auto res = PKDLog::range_report(tree, qs[i]);
                if (rep == 0) actual_cnt[i] = res.size();
            }
            ms += t.stop() * 1000.0;
        }
        std::cout << "[PkdLogTree]: range report time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        bool is_correct = true;
        for (size_t i = 0; i < qs.size(); i++) {
            if (actual_cnt[i] != cnt[i]) {
                std::cout << "[PKDLog ERROR] Query " << i << " failed. Expected: " << cnt[i] << ", Got: " << actual_cnt[i] << std::endl;
                is_correct = false; break;
            }
        }
        if (is_correct) std::cout << "[PkdLogTree] Accuracy: 100% (All " << qs.size() << " queries correct)" << std::endl;
    }
    
    inline void spatial_diff_test_latency(PT P, RQ querys, parlay::sequence<size_t>& batch_sizes, size_t insert_ratio, double p = 1.0) {
        auto P_conv = convert_points(P);
        auto tree0 = PKDLog::map_init(P_conv);
        auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());
        
        for (auto batch_size : batch_sizes) {
            std::cout << "[INFO] Batch Size: " << batch_size << std::endl;
            auto insert_num = batch_size / 10 * insert_ratio;
            auto delete_num = batch_size / 10 * (10 - insert_ratio);
            
            auto P_insert = P_insert_set.substr(0, insert_num);
            auto P_delete = P_delete_set.substr(0, delete_num);
            
            auto tree1 = PKDLog::map_delete(tree0, convert_points(P_delete), p);
            auto tree2 = PKDLog::map_insert(tree1, convert_points(P_insert), p);
            
            for (size_t i = 0; i < querys.size(); i++) {
                double avg_time = time_loop(
                    3, 1.0, 
                    [&](){},
                    [&](){
                        PKDLog::diff_type diff;
                        PKDLog::map_spatial_diff(tree0, tree2, querys[i], diff);
                    },
                    [&](){}
                );
                std::cout << std::fixed << std::setprecision(6) << i << " " << avg_time << std::endl;
            }
        }
    }
}