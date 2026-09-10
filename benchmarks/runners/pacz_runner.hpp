#pragma once
#include <string>
#include <vector>
#include "runner_common.hpp"
#include <silva/index/pacz.hpp>
#include <silva/core/global_config.hpp>

struct PACZRunner {
    static std::string name() { return "PACZ"; }
    
    PACZ::zmap master;
    std::vector<PACZ::zmap> history;
    
    PACZRunner() {}
    
    void build_base(const std::vector<Value>& P_base) {
        std::vector<geobase::Point> pts(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) pts[i] = geobase::Point(P_base[i].second, P_base[i].first.get<0>(), P_base[i].first.get<1>());
        
        // Dummy run
        std::vector<geobase::Point> sub_pts(pts.begin(), pts.begin() + std::min((size_t)10000, pts.size()));
        auto dummy = PACZ::map_init(sub_pts, false);
        
        master = PACZ::map_init(pts, false);
        history.push_back(master);
    }
    
    void reset_query_stats() {
        cpam::cpam_query_nodes_touched = 0;
    }
    
    size_t get_nodes_touched() const {
        return cpam::cpam_query_nodes_touched.load();
    }
    
    void range_query(geobase::Bounding_Box q_copy, parlay::sequence<geobase::Point>& shared_out, size_t& cnt, size_t& h) {
        cnt = PACZ::range_report(master, q_copy, shared_out, false);
        for(size_t j=0; j<cnt; j++) h += shared_out[j].id;
    }
    
    void knn_query(geobase::Point q, size_t k, size_t& cnt, size_t& h) {
        auto res = PACZ::knn(master, q, k);
        cnt = res.size();
        while(!res.empty()) { h += res.top().first.id; res.pop(); }
    }
    
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        auto v_pacz = PACZ::map_insert(adds, master);
        if (rems.size() > 0) v_pacz = PACZ::map_delete(rems, v_pacz);
        master = v_pacz;
        history.push_back(master);
    }
    
    std::pair<double, double> memory_usage() const {
        return {cpam::cpam_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0), 0.0};
    }
};
