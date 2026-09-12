#pragma once
#include <iostream>
#include <vector>
#include <silva/baselines/rlog_tree.hpp>
#include <helper/time_loop.h>

namespace RlogTest {
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
        auto run_f = [&]() {
            RlogTree tree(1);
            tree.build_base(P_conv);
            final_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
        };
        double ms = time_loop(5, 1.0, [](){}, run_f, [](){}) * 1000.0;
        std::cout << "[memory_MB]: " << final_mem << std::endl;
        std::cout << "[RlogTree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    void batch_insert_test(PT P_base, PT P_update, parlay::sequence<double>& batch_ratios, double p = 1.0) {
        auto P_conv = convert_points(P_base);
        auto n = P_base.size();
        auto rand_p = geobase::shuffle_point(P_update);
        parlay::parallel_for (0, rand_p.size(), [&](int i){
            rand_p[i].id = n + i;
        });
        auto P_update_conv = convert_points(rand_p);
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = P_update_conv.size() * ratio;
            if (cur_batch_size == 0) cur_batch_size = 1; // at least 1 point
            double ms = 0;
            double final_mem = 0;
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            for (int i = 0; i < 3; i++) {
                RlogTree tree(1);
                tree.build_base(P_conv);
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_update_conv.size()) {
                    size_t current_batch_size = std::min(cur_batch_size, P_update_conv.size() - offset);
                    std::vector<Value> batch(P_update_conv.begin() + offset, P_update_conv.begin() + offset + current_batch_size);
                    parlay::internal::timer batch_t;
                    tree.commit_inserts(batch);
                    tree.check_and_compact(p);
                    double b_time = batch_t.stop() * 1000.0;
                    if (i == 2) {
                        batch_times.push_back(b_time);
                        size_t rlog_bytes = (tree.insert_log.capacity() + tree.remove_log.capacity()) * sizeof(Value) + 
                                            tree.removed_ids.bucket_count() * sizeof(void*) + 
                                            tree.removed_ids.size() * sizeof(size_t);
                        batch_mems.push_back((boost_live_mem.load(std::memory_order_relaxed) + rlog_bytes) / (1024.0 * 1024.0));
                    }
                    offset += current_batch_size;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) {
                    size_t rlog_bytes = (tree.insert_log.capacity() + tree.remove_log.capacity()) * sizeof(Value) + 
                                        tree.removed_ids.bucket_count() * sizeof(void*) + 
                                        tree.removed_ids.size() * sizeof(size_t);
                    final_mem = (boost_live_mem.load(std::memory_order_relaxed) + rlog_bytes) / (1024.0 * 1024.0);
                }
            }
            std::cout << "[per_batch_time]: ";
            for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_batch_mem]: ";
            for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << final_mem << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
            std::cout << "[RlogTree]: batch insert time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void batch_delete_test(PT P_base, parlay::sequence<double>& batch_ratios, double p = 1.0) {
        auto rand_p = geobase::shuffle_point(P_base);
        auto P_conv = convert_points(P_base);
        auto P_delete_conv = convert_points(rand_p);
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = P_base.size() * ratio;
            if (cur_batch_size == 0) cur_batch_size = 1;
            double ms = 0;
            double final_mem = 0;
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            for (int i = 0; i < 3; i++) {
                RlogTree tree(1);
                tree.build_base(P_conv); // Build full tree first
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_delete_conv.size()) {
                    size_t current_batch_size = std::min(cur_batch_size, P_delete_conv.size() - offset);
                    std::vector<Value> batch(P_delete_conv.begin() + offset, P_delete_conv.begin() + offset + current_batch_size);
                    RlogBranch branch;
                    branch.remove_log = batch;
                    parlay::internal::timer batch_t;
                    tree.merge(branch); // Mark as deleted
                    tree.check_and_compact(p);
                    double b_time = batch_t.stop() * 1000.0;
                    if (i == 2) {
                        batch_times.push_back(b_time);
                        size_t rlog_bytes = (tree.insert_log.capacity() + tree.remove_log.capacity()) * sizeof(Value) + 
                                            tree.removed_ids.bucket_count() * sizeof(void*) + 
                                            tree.removed_ids.size() * sizeof(size_t);
                        batch_mems.push_back((boost_live_mem.load(std::memory_order_relaxed) + rlog_bytes) / (1024.0 * 1024.0));
                    }
                    offset += current_batch_size;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) {
                    size_t rlog_bytes = (tree.insert_log.capacity() + tree.remove_log.capacity()) * sizeof(Value) + 
                                        tree.removed_ids.bucket_count() * sizeof(void*) + 
                                        tree.removed_ids.size() * sizeof(size_t);
                    final_mem = (boost_live_mem.load(std::memory_order_relaxed) + rlog_bytes) / (1024.0 * 1024.0);
                }
            }
            std::cout << "[per_batch_time]: ";
            for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_batch_mem]: ";
            for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << final_mem << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
            std::cout << "[RlogTree]: batch delete time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void range_report_test(PT P, RQ querys, parlay::sequence<size_t>& cnt) {
        auto P_conv = convert_points(P);
        RlogTree tree(1);
        tree.build_base(P_conv);
        auto qs = querys;
        double ms = 0;
        
        std::vector<size_t> actual_cnt(qs.size(), 0);
        for (int rep = 0; rep < 3; rep++) {
            parlay::internal::timer t;
            for (size_t i = 0; i < qs.size(); i++) {
                auto res = tree.range_report(qs[i]);
                if (rep == 0) actual_cnt[i] = res.size();
            }
            ms += t.stop() * 1000.0;
        }
        std::cout << "[RlogTree]: range report time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        bool is_correct = true;
        for (size_t i = 0; i < qs.size(); i++) {
            if (actual_cnt[i] != cnt[i]) {
                std::cout << "[Rlog ERROR] Query " << i << " failed. Expected: " << cnt[i] << ", Got: " << actual_cnt[i] << std::endl;
                is_correct = false; break;
            }
        }
        if (is_correct) std::cout << "[RlogTree] Accuracy: 100% (All " << qs.size() << " queries correct)" << std::endl;
    }
}
