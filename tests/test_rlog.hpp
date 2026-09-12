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

    void batch_insert_test(PT P_base, PT P_update, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        auto P_update_conv = convert_points(P_update);
        for (double ratio : batch_ratios) {
            size_t chunk_size = P_update_conv.size() * ratio;
            if (chunk_size == 0) chunk_size = 1; // at least 1 point
            double ms = 0;
            double final_mem = 0;
            for (int i = 0; i < 3; i++) {
                RlogTree tree(1);
                tree.build_base(P_conv);
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_update_conv.size()) {
                    size_t current_chunk = std::min(chunk_size, P_update_conv.size() - offset);
                    std::vector<Value> batch(P_update_conv.begin() + offset, P_update_conv.begin() + offset + current_chunk);
                    tree.commit_inserts(batch);
                    offset += current_chunk;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) {
                    final_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
                }
            }
            std::cout << "[memory_MB]: " << final_mem << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
            std::cout << "[RlogTree]: batch insert time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void batch_delete_test(PT P_base, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        for (double ratio : batch_ratios) {
            size_t chunk_size = P_base.size() * ratio;
            if (chunk_size == 0) chunk_size = 1;
            double ms = 0;
            double final_mem = 0;
            for (int i = 0; i < 3; i++) {
                RlogTree tree(1);
                tree.build_base(P_conv); // Build full tree first
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_conv.size()) {
                    size_t current_chunk = std::min(chunk_size, P_conv.size() - offset);
                    std::vector<Value> batch(P_conv.begin() + offset, P_conv.begin() + offset + current_chunk);
                    RlogBranch branch;
                    branch.remove_log = batch;
                    tree.merge(branch); // Mark as deleted
                    offset += current_chunk;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) {
                    final_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
                }
            }
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
