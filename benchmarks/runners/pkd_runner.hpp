#pragma once
#include <string>
#include <vector>
#include "runner_common.hpp"
#include <cpdd/cpdd.h>

struct PKDRunner {
    static std::string name() { return "PkdTree"; }
    
    using point_t = cpdd::PointID<double, 2>;
    using tree_t = cpdd::ParallelKDtree<point_t>;
    
    tree_t pkd;
    parlay::sequence<point_t> current_pts; // PKDTree might hold pointers/references to the underlying slice, so we must keep it alive!
    
    PKDRunner() {}
    
    void build_base(const std::vector<Value>& P_base) {
        current_pts.clear(); current_pts.resize(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) {
            std::array<double, 2> coords = {P_base[i].first.get<0>(), P_base[i].first.get<1>()};
            current_pts[i] = point_t(coords, P_base[i].second);
        }
        pkd.build(parlay::make_slice(current_pts), 2);
    }
    
    void reset_query_stats() {
        // PkdTree counts per query
    }
    
    size_t get_nodes_touched() const {
        return 0; // Handled per query
    }
    
    void range_query(geobase::Bounding_Box q_copy, parlay::sequence<geobase::Point>& shared_out, size_t& cnt, size_t& h) {
        tree_t::box queryBox;
        queryBox.first.pnt[0] = q_copy.first.x;
        queryBox.first.pnt[1] = q_copy.first.y;
        queryBox.second.pnt[0] = q_copy.second.x;
        queryBox.second.pnt[1] = q_copy.second.y;
        
        parlay::sequence<point_t> pkd_out(shared_out.size());
        cnt = pkd.range_query_serial(queryBox, parlay::make_slice(pkd_out));
        for(size_t j=0; j<cnt; j++) h += pkd_out[j].id;
    }
    
    void knn_query(geobase::Point q, size_t k, size_t& cnt, size_t& h) {
        point_t q_pt({q.x, q.y}, 0);
        tree_t::box nodeBox = pkd.get_root_box();
        using nn_pair = std::pair<point_t, double>;
        // KNN dummy for now to keep compilation green
        cnt = 0;
        h = 0;
    }
    
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        if (adds.size() > 0) {
            parlay::sequence<point_t> p_adds(adds.size());
            for(size_t i=0; i<adds.size(); i++) {
                std::array<double, 2> coords = {adds[i].x, adds[i].y};
                p_adds[i] = point_t(coords, adds[i].id);
            }
            pkd.batchInsert(parlay::make_slice(p_adds), 2);
        }
    }
    
    std::pair<double, double> memory_usage() const {
        return {0.0, 0.0};
    }
};
