#pragma once
#include <iostream>
#include <vector>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <psi/base_tree.h>
#include <psi/dependence/splitter.h>
#include <psi/p_tree.h>
#include "test_utils.hpp"

namespace SPACTest {

struct spac_aug_id_code {
	using id_type = int_fast32_t;
	using curve_code_type = uint64_t;
	curve_code_type code;
	id_type id;
	spac_aug_id_code() : code(0), id(0) {}
	void set_member(curve_code_type const &val) { code = val; }
	bool operator<(spac_aug_id_code const &rhs) const {
		return code == rhs.code ? id < rhs.id : code < rhs.code;
	}
	bool operator==(spac_aug_id_code const &rhs) const { return id == rhs.id; }
};

using Point = psi::aug_point<double, 2, spac_aug_id_code>;
using points_type = parlay::sequence<Point>;
using SplitRule = psi::spatial_filling_curve<psi::morton_curve<Point>>;
using Tree = psi::p_tree<psi::tree_traits<Point, SplitRule>>;

inline Tree map_init(parlay::sequence<geobase::Point>& P_base) {
    auto n = P_base.size();
    points_type pts(n);
    parlay::parallel_for(0, n, [&](int i){
        pts[i][0] = P_base[i].x;
        pts[i][1] = P_base[i].y;
        pts[i].aug.id = P_base[i].id;
    });
    Tree tree;
    tree.build(std::move(pts));
    return tree;
}

inline void batch_insert_test(parlay::sequence<geobase::Point>& P_base, parlay::sequence<geobase::Point>& P_update, parlay::sequence<double> &batch_ratios, bool single_version = true) {
    auto m1 = map_init(P_base);
    auto rand_p = geobase::shuffle_point(P_update);
    auto n = P_base.size();
    parlay::parallel_for(0, rand_p.size(), [&](int i){
        rand_p[i].id = n + i;
    });

    for (auto ratio: batch_ratios){
        size_t cur_batch_size = std::max<size_t>(1, rand_p.size() * ratio);
        std::cout << "[Testing Ratio]: " << ratio << std::endl;
        
        std::vector<double> batch_times;
        std::vector<double> batch_mems;
        double total_ms = 0;

        Tree sv_current_tree;
        double total_avg = time_loop(
            3, 1.0,
            [&]() {
                if (single_version) {
                    sv_current_tree = map_init(P_base);
                } else {
                    std::cerr << "SPaCtree does not support multi-version!" << std::endl;
                    exit(1);
                }
                batch_times.clear();
                batch_mems.clear();
            },
            [&]() {
                for (size_t i = 0; i < rand_p.size(); i += cur_batch_size) {
                    size_t current_batch_size = std::min(cur_batch_size, rand_p.size() - i);
                    auto P2 = rand_p.substr(i, current_batch_size);

                    points_type adds(current_batch_size);
                    parlay::parallel_for(0, current_batch_size, [&](int j){
                        adds[j][0] = P2[j].x;
                        adds[j][1] = P2[j].y;
                        adds[j].aug.id = P2[j].id;
                    });

                    parlay::internal::timer t;
                    sv_current_tree.batch_insert(parlay::make_slice(adds));
                    batch_times.push_back(t.next_time() * 1000.0);
                    
                    double mem_mb = 0;
                    batch_mems.push_back(mem_mb);
                }
            },
            [&]() {
                sv_current_tree.delete_tree();
            }
        );

        total_ms = total_avg * 1000.0;

        std::cout << "[per_batch_time]: ";
        for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
        std::cout << std::endl;
        std::cout << "[per_batch_mem]: ";
        for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
        std::cout << std::endl;
        std::cout << "[memory_MB]: " << batch_mems.back() << std::endl;
        std::cout << "[batch_ratio]: " << ratio << std::endl;
        std::cout << "[SPaC-tree]: batch insert time (avg): " << (total_ms / 1000.0) << std::endl;
    }
}

inline void batch_delete_test(parlay::sequence<geobase::Point>& P_base, parlay::sequence<double> &batch_ratios, bool single_version = true) {
    auto rand_p = geobase::shuffle_point(P_base);

    for (auto ratio: batch_ratios){
        size_t cur_batch_size = std::max<size_t>(1, rand_p.size() * ratio);
        std::cout << "[Testing Ratio]: " << ratio << std::endl;
        
        std::vector<double> batch_times;
        std::vector<double> batch_mems;
        double total_ms = 0;

        Tree sv_current_tree;
        double total_avg = time_loop(
            3, 1.0,
            [&]() {
                if (single_version) {
                    sv_current_tree = map_init(P_base);
                } else {
                    std::cerr << "SPaCtree does not support multi-version!" << std::endl;
                    exit(1);
                }
                batch_times.clear();
                batch_mems.clear();
            },
            [&]() {
                for (size_t i = 0; i < rand_p.size(); i += cur_batch_size) {
                    size_t current_batch_size = std::min(cur_batch_size, rand_p.size() - i);
                    auto P2 = rand_p.substr(i, current_batch_size);

                    points_type adds(current_batch_size);
                    parlay::parallel_for(0, current_batch_size, [&](int j){
                        adds[j][0] = P2[j].x;
                        adds[j][1] = P2[j].y;
                        adds[j].aug.id = P2[j].id;
                    });

                    parlay::internal::timer t;
                    sv_current_tree.batch_delete(parlay::make_slice(adds));
                    batch_times.push_back(t.next_time() * 1000.0);
                    
                    double mem_mb = 0;
                    batch_mems.push_back(mem_mb);
                }
            },
            [&]() {
                sv_current_tree.delete_tree();
            }
        );

        total_ms = total_avg * 1000.0;

        std::cout << "[per_batch_time]: ";
        for(size_t j = 0; j < batch_times.size(); j++) std::cout << batch_times[j] << (j==batch_times.size()-1 ? "" : ",");
        std::cout << std::endl;
        std::cout << "[per_batch_mem]: ";
        for(size_t j = 0; j < batch_mems.size(); j++) std::cout << batch_mems[j] << (j==batch_mems.size()-1 ? "" : ",");
        std::cout << std::endl;
        std::cout << "[memory_MB]: " << batch_mems.back() << std::endl;
        std::cout << "[batch_ratio]: " << ratio << std::endl;
        std::cout << "[SPaC-tree]: batch delete time (avg): " << (total_ms / 1000.0) << std::endl;
    }
}

}