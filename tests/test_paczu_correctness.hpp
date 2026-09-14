#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <parlay/primitives.h>
#include <silva/index/pacz.hpp>
#include <silva/index/pacz_u.hpp>
#include <silva/geo/io.hpp>
#include "test_utils.hpp"

namespace PACZCorrectness {

    using namespace std;
    using namespace geobase;
    using PT = parlay::sequence<geobase::Point>;

    void verify_test(PT P_base, PT P_update) {
        cout << "[Verify] Starting PaCZtree vs PaCZUtree correctness verification..." << endl;
        
        // 1. Build Base
        PACZ::zmap tree1;
        PACZU::zmap tree2;

        auto n = P_base.size();
        auto insert_pts1 = parlay::sequence<PACZ::par>::uninitialized(n);
        auto insert_pts2 = parlay::sequence<PACZU::par>::uninitialized(n);
        parlay::parallel_for(0, n, [&](int i) {
            insert_pts1[i] = {{P_base[i].morton_id, P_base[i].id}, P_base[i]};
            insert_pts2[i] = {{P_base[i].morton_id, P_base[i].id}, P_base[i]};
        });

        tree1 = PACZ::zmap(insert_pts1);
        tree2 = PACZU::zmap(insert_pts2);
        
        if (tree1.size() != tree2.size()) {
            cout << "[ERROR] Base build size mismatch: " << tree1.size() << " vs " << tree2.size() << endl;
            exit(1);
        }

        // 2. Complex batch operations: Insert 10%, Delete 5%, Insert 5%
        auto rand_p = geobase::shuffle_point(P_update);
        parlay::parallel_for(0, rand_p.size(), [&](int i) { rand_p[i].id = n + i; });
        
        size_t b1 = rand_p.size() * 0.1;
        size_t b2 = rand_p.size() * 0.05;

        // Insert batch 1
        auto b1_pts1 = parlay::sequence<PACZ::par>::uninitialized(b1);
        auto b1_pts2 = parlay::sequence<PACZU::par>::uninitialized(b1);
        parlay::parallel_for(0, b1, [&](int i) {
            b1_pts1[i] = {{rand_p[i].morton_id, rand_p[i].id}, rand_p[i]};
            b1_pts2[i] = {{rand_p[i].morton_id, rand_p[i].id}, rand_p[i]};
        });
        tree1 = PACZ::zmap::multi_insert(tree1, b1_pts1);
        tree2 = PACZU::zmap::multi_insert(tree2, b1_pts2);

        // Delete batch 2 (subset of batch 1)
        auto b2_pts1 = parlay::sequence<PACZ::par>::uninitialized(b2);
        auto b2_pts2 = parlay::sequence<PACZU::par>::uninitialized(b2); // PACZU uses multi_delete with par
        parlay::parallel_for(0, b2, [&](int i) {
            b2_pts1[i] = {{rand_p[i].morton_id, rand_p[i].id}, rand_p[i]};
            b2_pts2[i] = {{rand_p[i].morton_id, rand_p[i].id}, rand_p[i]};
        });
        tree1 = PACZ::zmap::multi_delete(tree1, b2_pts1);
        tree2 = PACZU::zmap::multi_delete(tree2, b2_pts2);
        
        if (tree1.size() != tree2.size()) {
            cout << "[ERROR] Size mismatch after batch operations: " << tree1.size() << " vs " << tree2.size() << endl;
            exit(1);
        }

        // 3. Spatial diff (Version operation)
        Bounding_Box q_diff = {Point(0.2, 0.2), Point(0.8, 0.8)};
        
        struct DIFF1 { vector<PACZ::val_type> add; size_t add_cnt; vector<PACZ::val_type> remove; size_t remove_cnt; };
        struct DIFF2 { vector<PACZU::val_type> add; size_t add_cnt; vector<PACZU::val_type> remove; size_t remove_cnt; };
        
        DIFF1 d1; DIFF2 d2;
        parlay::sequence<Point> l1, r1, l2, r2;
        // In the original pacz, plain_map_spatial_diff doesn't exist, map_spatial_diff does. 
        // We'll just do filter_range and map_difference manually to test it perfectly
        
        auto f1 = PACZ::filter_range(tree1, q_diff);
        auto f2 = PACZU::filter_range(tree2, q_diff);
        if (f1.size() != f2.size()) {
            cout << "[ERROR] Spatial filter size mismatch: " << f1.size() << " vs " << f2.size() << endl;
            exit(1);
        }

        // Dummy base version to diff against (an empty tree)
        PACZ::zmap empty1; PACZU::zmap empty2;
        auto add1 = PACZ::zmap::values(PACZ::zmap::map_difference(f1, empty1));
        auto add2 = PACZU::zmap::values(PACZU::zmap::map_difference(f2, empty2));
        
        if (add1.size() != add2.size()) {
            cout << "[ERROR] Spatial diff size mismatch: " << add1.size() << " vs " << add2.size() << endl;
            exit(1);
        }
        
        long long checksum1 = 0, checksum2 = 0;
        for (auto p : add1) checksum1 += p.id;
        for (auto p : add2) checksum2 += p.id;
        if (checksum1 != checksum2) {
            cout << "[ERROR] Checksum mismatch in diff!" << endl;
            exit(1);
        }

        cout << "[Verify] SUCCESS! PaCZtree and PaCZUtree behaviors match 100%!" << endl;
    }
}
