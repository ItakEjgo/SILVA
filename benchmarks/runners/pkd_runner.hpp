#pragma once
#include <string>
#include <vector>
#include "runner_common.hpp"
#include <psi/kd_tree.h>
#include <psi/dependence/splitter.h>

struct pkd_aug_id {
    using id_type = int;
    id_type id;

    bool operator<(pkd_aug_id const &rhs) const { return id < rhs.id; }
    bool operator==(pkd_aug_id const &rhs) const { return id == rhs.id; }
};

struct PKDRunner {
    static std::string name() { return "PkdTree"; }
    
    using point_t = psi::aug_point<double, 2, pkd_aug_id>;
    using SplitRule = psi::orthogonal_split_rule<psi::max_stretch_dim<point_t>, psi::object_median<point_t>>;
    using tree_t = psi::kd_tree<psi::tree_traits<point_t, SplitRule>>;
    
    tree_t pkd;
    parlay::sequence<point_t> current_pts; // PKDTree might hold pointers/references to the underlying slice, so we must keep it alive!
    
    PKDRunner() {}
    ~PKDRunner() { pkd.delete_tree(); }
    
    void build_base(const std::vector<Value>& P_base) {
        current_pts.clear(); current_pts.resize(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) {
            current_pts[i][0] = P_base[i].first.get<0>();
            current_pts[i][1] = P_base[i].first.get<1>();
            current_pts[i].aug.id = P_base[i].second;
        }
        pkd.build(parlay::make_slice(current_pts));
    }
    
    void reset_query_stats() {
        // PkdTree counts per query
    }
    
    size_t get_nodes_touched() const {
        return 0; // Handled per query
    }
    
    void range_query(geobase::Bounding_Box q_copy, parlay::sequence<geobase::Point>& shared_out, size_t& cnt, size_t& h) {
        tree_t::box_type queryBox;
        queryBox.first[0] = q_copy.first.x;
        queryBox.first[1] = q_copy.first.y;
        queryBox.second[0] = q_copy.second.x;
        queryBox.second[1] = q_copy.second.y;
        
        parlay::sequence<point_t> pkd_out(shared_out.size());
        auto res = pkd.range_query(queryBox, parlay::make_slice(pkd_out));
        cnt = res.first;
        for(size_t j=0; j<cnt; j++) h += pkd_out[j].aug.id;
    }
    
    void knn_query(geobase::Point q, size_t k, size_t& cnt, size_t& h) {
        point_t q_pt;
        q_pt[0] = q.x;
        q_pt[1] = q.y;
        q_pt.aug.id = -1;

        using dis_type = typename point_t::dis_type;
        using nn_pair = std::pair<point_t, dis_type>;

        parlay::sequence<nn_pair> knn_result(k, nn_pair(q_pt, 0));
        psi::bounded_queue<point_t, nn_pair> bq(parlay::make_slice(knn_result));
        
        pkd.knn(q_pt, bq);
        
        cnt = k; 
        for (size_t i = 0; i < k; i++) h += knn_result[i].first.aug.id;
    }
    
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        if (adds.size() > 0) {
            parlay::sequence<point_t> p_adds(adds.size());
            for(size_t i=0; i<adds.size(); i++) {
                p_adds[i][0] = adds[i].x;
                p_adds[i][1] = adds[i].y;
                p_adds[i].aug.id = adds[i].id;
            }
            pkd.batch_insert(parlay::make_slice(p_adds));
        }
        if (rems.size() > 0) {
            parlay::sequence<point_t> p_rems(rems.size());
            for(size_t i=0; i<rems.size(); i++) {
                p_rems[i][0] = rems[i].x;
                p_rems[i][1] = rems[i].y;
                p_rems[i].aug.id = rems[i].id;
            }
            pkd.batch_delete(parlay::make_slice(p_rems));
        }
    }
    
    size_t calculate_tree_memory(auto T) const {
        if (T == nullptr) return 0;
        if (T->is_leaf) {
            return sizeof(typename tree_t::leaf_type) + tree_t::traits_type::leaf_capacity * sizeof(point_t);
        }
        
        auto TI = static_cast<typename tree_t::interior_type*>(T);
        size_t l = 0, r = 0;
        parlay::par_do_if(T->size > 1000,
            [&]() { l = calculate_tree_memory(TI->left); },
            [&]() { r = calculate_tree_memory(TI->right); }
        );
        return sizeof(typename tree_t::interior_type) + l + r;
    }

    std::pair<double, double> memory_usage() const {
        double mem_mb = calculate_tree_memory(const_cast<tree_t&>(pkd).get_root()) / (1024.0 * 1024.0);
        return {mem_mb, mem_mb};
    }
};
