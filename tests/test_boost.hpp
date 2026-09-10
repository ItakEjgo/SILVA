#pragma once
#include <iostream>
#include <vector>
#include "../benchmarks/runners/boost_runner.hpp"
#include <helper/time_loop.h>

namespace BoostTest {
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
            BoostRunner runner;
            runner.build_base(P_conv);
        };
        double ms = time_loop(5, 1.0, [](){}, run_f, [](){}) * 1000.0;
        std::cout << "[BoostRtree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    void batch_insert_test(PT P, parlay::sequence<size_t>& batch_sizes) {
        auto P_conv = convert_points(P);
        for (size_t b_size : batch_sizes) {
            size_t actual_b = std::min(b_size, P_conv.size());
            parlay::sequence<geobase::Point> adds(actual_b);
            for(size_t i=0; i<actual_b; i++) adds[i] = P[i];
            parlay::sequence<geobase::Point> rems; // empty

            double ms = 0;
            for (int i = 0; i < 3; i++) {
                BoostRunner runner;
                parlay::internal::timer t;
                runner.commit(adds, rems);
                ms += t.stop() * 1000.0;
            }
            std::cout << "[batch_size]: " << actual_b << std::endl;
            std::cout << "[BoostRtree]: batch insert time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void batch_delete_test(PT P, parlay::sequence<size_t>& batch_sizes) {
        auto P_conv = convert_points(P);
        for (size_t b_size : batch_sizes) {
            size_t actual_b = std::min(b_size, P_conv.size());
            parlay::sequence<geobase::Point> adds; // empty
            parlay::sequence<geobase::Point> rems(actual_b);
            for(size_t i=0; i<actual_b; i++) rems[i] = P[i];

            double ms = 0;
            for (int i = 0; i < 3; i++) {
                BoostRunner runner;
                runner.build_base(P_conv);
                parlay::internal::timer t;
                runner.commit(adds, rems);
                ms += t.stop() * 1000.0;
            }
            std::cout << "[batch_size]: " << actual_b << std::endl;
            std::cout << "[BoostRtree]: batch delete time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void range_report_test(PT P, RQ qs, parlay::sequence<size_t>& cnt) {
        auto P_conv = convert_points(P);
        BoostRunner runner;
        runner.build_base(P_conv);
        
        parlay::sequence<geobase::Point> shared_out(P.size());
        std::vector<size_t> actual_cnt(qs.size(), 0);
        
        double ms = 0;
        for (int rep = 0; rep < 3; rep++) {
            parlay::internal::timer t;
            for (size_t i = 0; i < qs.size(); i++) {
                size_t h = 0;
                runner.range_query(qs[i], shared_out, actual_cnt[i], h);
            }
            ms += t.stop() * 1000.0;
        }

        std::cout << "[BoostRtree]: range report time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        
        bool is_correct = true;
        for (size_t i = 0; i < qs.size(); i++) {
            if (actual_cnt[i] != cnt[i]) {
                std::cout << "[Boost ERROR] Query " << i << " failed. Expected: " << cnt[i] << ", Got: " << actual_cnt[i] << std::endl;
                is_correct = false; break;
            }
        }
        if (is_correct) std::cout << "[BoostRtree] Accuracy: 100% (All " << qs.size() << " queries correct)" << std::endl;
    }
}
