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
        auto run_f = [&]() {
            RlogTree tree(1);
            tree.build_base(P_conv);
        };
        double ms = time_loop(5, 1.0, [](){}, run_f, [](){}) * 1000.0;
        std::cout << "[RlogTree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    void batch_insert_test(PT P, parlay::sequence<size_t>& batch_sizes) {
        auto P_conv = convert_points(P);
        for (size_t b_size : batch_sizes) {
            size_t actual_b = std::min(b_size, P_conv.size());
            std::vector<Value> batch(P_conv.begin(), P_conv.begin() + actual_b);
            double ms = 0;
            for (int i = 0; i < 3; i++) {
                RlogTree tree(1);
                parlay::internal::timer t;
                tree.commit_inserts(batch);
                ms += t.stop() * 1000.0;
            }
            std::cout << "[batch_size]: " << actual_b << std::endl;
            std::cout << "[RlogTree]: batch insert time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void batch_delete_test(PT P, parlay::sequence<size_t>& batch_sizes) {
        auto P_conv = convert_points(P);
        for (size_t b_size : batch_sizes) {
            size_t actual_b = std::min(b_size, P_conv.size());
            std::vector<Value> batch(P_conv.begin(), P_conv.begin() + actual_b);
            double ms = 0;
            for (int i = 0; i < 3; i++) {
                RlogTree tree(1);
                tree.build_base(P_conv); // Build full tree first
                parlay::internal::timer t;
                RlogBranch branch;
                branch.remove_log = batch;
                tree.merge(branch); // Mark as deleted
                ms += t.stop() * 1000.0;
            }
            std::cout << "[batch_size]: " << actual_b << std::endl;
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
