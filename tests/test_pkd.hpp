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
        auto run_f = [&]() {
            PKDRunner runner;
            runner.build_base(P_conv);
        };
        double ms = time_loop(5, 1.0, [](){}, run_f, [](){}) * 1000.0;
        std::cout << "[PkdTree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    void batch_insert_test(PT P_base, PT P_update, parlay::sequence<size_t>& batch_sizes) {
        auto P_conv = convert_points(P_base);
        for (size_t b_size : batch_sizes) {
            size_t actual_b = std::min(b_size, P_update.size());
            parlay::sequence<geobase::Point> adds(actual_b);
            for(size_t i=0; i<actual_b; i++) adds[i] = P_update[i];
            parlay::sequence<geobase::Point> rems; // empty

            double ms = 0;
            for (int i = 0; i < 3; i++) {
                PKDRunner runner;
                runner.build_base(P_conv);
                parlay::internal::timer t;
                runner.commit(adds, rems);
                ms += t.stop() * 1000.0;
            }
            std::cout << "[batch_size]: " << actual_b << std::endl;
            std::cout << "[PkdTree]: batch insert time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void batch_delete_test(PT P_base, parlay::sequence<size_t>& batch_sizes) {
        auto P_conv = convert_points(P_base);
        for (size_t b_size : batch_sizes) {
            size_t actual_b = std::min(b_size, P_base.size());
            parlay::sequence<geobase::Point> adds; // empty
            parlay::sequence<geobase::Point> rems(actual_b);
            for(size_t i=0; i<actual_b; i++) rems[i] = P_base[i];

            double ms = 0;
            for (int i = 0; i < 3; i++) {
                PKDRunner runner;
                runner.build_base(P_conv);
                parlay::internal::timer t;
                runner.commit(adds, rems);
                ms += t.stop() * 1000.0;
            }
            std::cout << "[batch_size]: " << actual_b << std::endl;
            std::cout << "[PkdTree]: batch delete time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
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
