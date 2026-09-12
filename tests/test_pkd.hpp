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
        auto run_f = [&]() {
            PKDRunner runner;
            runner.build_base(P_conv);
            final_mem = runner.memory_usage().first;
        };
        double ms = time_loop(5, 1.0, [](){}, run_f, [](){}) * 1000.0;
        std::cout << "[memory_MB]: " << final_mem << std::endl;
        std::cout << "[PkdTree]: build time (avg): " << ms / 1000.0 << std::endl;
    }

    void batch_insert_test(PT P_base, PT P_update, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        for (double ratio : batch_ratios) {
            size_t chunk_size = P_update.size() * ratio;
            if (chunk_size == 0) chunk_size = 1;
            
            double ms = 0;
            double final_mem = 0;
            std::vector<double> chunk_times;
            std::vector<double> chunk_mems;
            for (int i = 0; i < 3; i++) {
                PKDRunner runner;
                runner.build_base(P_conv);
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_update.size()) {
                    size_t current_chunk = std::min(chunk_size, P_update.size() - offset);
                    parlay::sequence<geobase::Point> adds(current_chunk);
                    for(size_t j=0; j<current_chunk; j++) adds[j] = P_update[offset + j];
                    parlay::sequence<geobase::Point> rems; // empty
                    parlay::internal::timer chunk_t;

                    runner.commit(adds, rems);

                    double c_time = chunk_t.stop() * 1000.0;

                    if (i == 2) {

                        chunk_times.push_back(c_time);

                        chunk_mems.push_back(runner.memory_usage().first);

                    }

                    offset += current_chunk;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) final_mem = runner.memory_usage().first;
            }
            std::cout << "[per_chunk_time]: ";
            for(size_t j = 0; j < chunk_times.size(); j++) std::cout << chunk_times[j] << (j==chunk_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_chunk_mem]: ";
            for(size_t j = 0; j < chunk_mems.size(); j++) std::cout << chunk_mems[j] << (j==chunk_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << final_mem << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
            std::cout << "[PkdTree]: batch insert time (avg): " << (ms / 3.0) / 1000.0 << std::endl;
        }
    }

    void batch_delete_test(PT P_base, parlay::sequence<double>& batch_ratios) {
        auto P_conv = convert_points(P_base);
        for (double ratio : batch_ratios) {
            size_t chunk_size = P_base.size() * ratio;
            if (chunk_size == 0) chunk_size = 1;
            
            double ms = 0;
            double final_mem = 0;
            std::vector<double> chunk_times;
            std::vector<double> chunk_mems;
            for (int i = 0; i < 3; i++) {
                PKDRunner runner;
                runner.build_base(P_conv);
                
                parlay::internal::timer t;
                size_t offset = 0;
                while (offset < P_base.size()) {
                    size_t current_chunk = std::min(chunk_size, P_base.size() - offset);
                    parlay::sequence<geobase::Point> adds; // empty
                    parlay::sequence<geobase::Point> rems(current_chunk);
                    for(size_t j=0; j<current_chunk; j++) rems[j] = P_base[offset + j];
                    parlay::internal::timer chunk_t;

                    runner.commit(adds, rems);

                    double c_time = chunk_t.stop() * 1000.0;

                    if (i == 2) {

                        chunk_times.push_back(c_time);

                        chunk_mems.push_back(runner.memory_usage().first);

                    }

                    offset += current_chunk;
                }
                ms += t.stop() * 1000.0;
                if (i == 2) final_mem = runner.memory_usage().first;
            }
            std::cout << "[per_chunk_time]: ";
            for(size_t j = 0; j < chunk_times.size(); j++) std::cout << chunk_times[j] << (j==chunk_times.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[per_chunk_mem]: ";
            for(size_t j = 0; j < chunk_mems.size(); j++) std::cout << chunk_mems[j] << (j==chunk_mems.size()-1 ? "" : ",");
            std::cout << std::endl;
            std::cout << "[memory_MB]: " << final_mem << std::endl;
            std::cout << "[batch_ratio]: " << ratio << std::endl;
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
