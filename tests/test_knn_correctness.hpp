#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include "../benchmarks/runners/mvq_runner.hpp"
#include "../benchmarks/runners/pacz_runner.hpp"
#include "../benchmarks/runners/boost_runner.hpp"
#include "../benchmarks/runners/rlog_runner.hpp"
#include "../benchmarks/runners/pkd_runner.hpp"
#include "../benchmarks/runners/pkd_log_runner.hpp"

namespace KNNVerify {
    using PT = parlay::sequence<geobase::Point>;

    void run_verification(PT P) {
        std::vector<Value> P_conv(P.size());
        for (size_t i = 0; i < P.size(); i++) {
            P_conv[i] = std::make_pair(BoostPoint(P[i].x, P[i].y), P[i].id);
        }

        std::cout << "--- Building Trees for Validation ---" << std::endl;
        MVQRunner mvq_r(mvq::Config::get().leaf_size); mvq_r.build_base(P_conv);
        PACZRunner pacz_r; pacz_r.build_base(P_conv);
        BoostRunner boost_r; boost_r.build_base(P_conv);
        RlogRunner rlog_r(1); rlog_r.build_base(P_conv);
        PKDRunner pkd_r; pkd_r.build_base(P_conv);
        PKDLogRunner pkd_log_r; pkd_log_r.build_base(P_conv);

        // Generate 100 random queries
        srand(42);
        std::vector<geobase::Point> queries(100);
        for(int i=0; i<100; i++) {
            queries[i] = P[rand() % P.size()];
        }

        std::vector<int> Ks = {1, 10, 100};
        
        for (int K : Ks) {
            std::cout << "\n=== Running " << K << "NN Correctness Check (100 Queries) ===" << std::endl;
            size_t mvq_hash = 0, pacz_hash = 0, boost_hash = 0, rlog_hash = 0, pkd_hash = 0, pkd_log_hash = 0;
            

            double mvq_time = 0, pacz_time = 0, boost_time = 0, rlog_time = 0, pkd_time = 0, pkd_log_time = 0;
            
            for (auto& q : queries) {
                size_t cnt = 0, h = 0;
                parlay::internal::timer t;
                
                t.start(); mvq_r.knn_query(q, K, cnt, h); mvq_time += t.stop(); mvq_hash += h;
                
                h = 0; t.start(); pacz_r.knn_query(q, K, cnt, h); pacz_time += t.stop(); pacz_hash += h;
                
                h = 0; t.start(); boost_r.knn_query(q, K, cnt, h); boost_time += t.stop(); boost_hash += h;
                
                h = 0; t.start(); rlog_r.knn_query(q, K, cnt, h); rlog_time += t.stop(); rlog_hash += h;
                
                h = 0; t.start(); pkd_r.knn_query(q, K, cnt, h); pkd_time += t.stop(); pkd_hash += h;
                
                h = 0; t.start(); pkd_log_r.knn_query(q, K, cnt, h); pkd_log_time += t.stop(); pkd_log_hash += h;
            }
            
            std::cout << std::left << std::setw(15) << "Algorithm" << std::setw(25) << "Checksum (Sum of IDs)" << "Time (ms)" << std::endl;
            std::cout << "--------------------------------------------------------" << std::endl;
            std::cout << std::left << std::setw(15) << "MVQ" << std::setw(25) << mvq_hash << (mvq_time * 1000.0) << std::endl;
            std::cout << std::left << std::setw(15) << "PACZ" << std::setw(25) << pacz_hash << (pacz_time * 1000.0) << std::endl;
            std::cout << std::left << std::setw(15) << "BoostRTree" << std::setw(25) << boost_hash << (boost_time * 1000.0) << std::endl;
            std::cout << std::left << std::setw(15) << "RlogTree" << std::setw(25) << rlog_hash << (rlog_time * 1000.0) << std::endl;
            std::cout << std::left << std::setw(15) << "PkdTree" << std::setw(25) << pkd_hash << (pkd_time * 1000.0) << std::endl;
            std::cout << std::left << std::setw(15) << "PkdLogTree" << std::setw(25) << pkd_log_hash << (pkd_log_time * 1000.0) << std::endl;
if (mvq_hash == pacz_hash && pacz_hash == boost_hash && boost_hash == rlog_hash && rlog_hash == pkd_hash && pkd_hash == pkd_log_hash) {
                std::cout << ">> RESULT: [PASS] All 6 algorithms returned identical results!" << std::endl;
            } else {
                std::cout << ">> RESULT: [FAIL] Checksums mismatch! Correctness bug detected." << std::endl;
            }
        }
    }
}
