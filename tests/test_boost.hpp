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
        double final_mem = 0;
        auto run_f = [&]() {
            BoostRunner runner;
            runner.build_base(P_conv);
            final_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
        };
        double ms = time_loop(5, 1.0, [](){}, run_f, [](){}) * 1000.0;
        std::cout << "[memory_MB]: " << final_mem << std::endl;
        std::cout << "[BoostRtree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    void batch_insert_test(PT P_base, PT P_update, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = P_update.size() * ratio;
            if (cur_batch_size == 0) cur_batch_size = 1;
            
            double ms = 0;
            double final_mem = 0;
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            for (int i = 0; i < 3; i++) {
                BoostRunner runner;
                runner.build_base(P_conv);
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_update.size()) {
                    size_t current_batch_size = std::min(cur_batch_size, P_update.size() - offset);
                    parlay::sequence<geobase::Point> adds(current_batch_size);
                    for(size_t j=0; j<current_batch_size; j++) adds[j] = P_update[offset + j];
                    parlay::sequence<geobase::Point> rems; // empty
                    parlay::internal::timer batch_t;

                    runner.commit(adds, rems);

                    double b_time = batch_t.stop() * 1000.0;

                    if (i == 2) {

                        batch_times.push_back(b_time);

                        batch_mems.push_back(runner.memory_usage().first);

                    }

                    offset += current_batch_size;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) {
                    final_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
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
            std::cout << "[BoostRtree]: batch insert time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void batch_delete_test(PT P_base, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        for (double ratio : batch_ratios) {
            size_t cur_batch_size = P_base.size() * ratio;
            if (cur_batch_size == 0) cur_batch_size = 1;
            
            double ms = 0;
            double final_mem = 0;
            std::vector<double> batch_times;
            std::vector<double> batch_mems;
            for (int i = 0; i < 3; i++) {
                BoostRunner runner;
                runner.build_base(P_conv);
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_base.size()) {
                    size_t current_batch_size = std::min(cur_batch_size, P_base.size() - offset);
                    parlay::sequence<geobase::Point> adds; // empty
                    parlay::sequence<geobase::Point> rems(current_batch_size);
                    for(size_t j=0; j<current_batch_size; j++) rems[j] = P_base[offset + j];
                    parlay::internal::timer batch_t;

                    runner.commit(adds, rems);

                    double b_time = batch_t.stop() * 1000.0;

                    if (i == 2) {

                        batch_times.push_back(b_time);

                        batch_mems.push_back(runner.memory_usage().first);

                    }

                    offset += current_batch_size;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) {
                    final_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
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
