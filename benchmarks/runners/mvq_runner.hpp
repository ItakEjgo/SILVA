#pragma once
#include <string>
#include <vector>
#include "runner_common.hpp"
#include <silva/index/mvq.hpp>
#include <silva/core/global_config.hpp>

struct MVQRunner {
    static std::string name() { return "MVQ"; }
    
    mvq::Tree* tree;
    std::shared_ptr<mvq::BaseNode> master;
    std::vector<std::shared_ptr<mvq::BaseNode>> history;
    
    MVQRunner(int leaf_size) {
        tree = new mvq::Tree(leaf_size);
    }
    
    ~MVQRunner() {
        delete tree;
    }
    
    void build_base(const std::vector<Value>& P_base) {
        // MVQ takes vector<Point> directly from P_base_set
        std::vector<geobase::Point> pts(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) pts[i] = geobase::Point(P_base[i].first.get<0>(), P_base[i].first.get<1>(), P_base[i].second);
        tree->build(pts);
        master = tree->root;
        history.push_back(master);
    }
    
    void reset_query_stats() {
        mvq::query_nodes_touched = 0;
    }
    
    size_t get_nodes_touched() const {
        return mvq::query_nodes_touched.load();
    }
    
    void range_query(geobase::Bounding_Box q_copy, parlay::sequence<geobase::Point>& shared_out, size_t& cnt, size_t& h) {
        cnt = 0;
        tree->range_report(master, q_copy, mvq::Config::get().largest_mbr, cnt, shared_out);
        for(size_t j=0; j<cnt; j++) h += shared_out[j].id;
    }
    void knn_query(const geobase::Point& q, size_t k, size_t& cnt, size_t& h) {
        auto res = tree->knn_report(k, q, mvq::Config::get().largest_mbr);
        cnt = res.size();
        while(!res.empty()) { h += res.top().first.id; res.pop(); }
    }
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        master = tree->commit(master, adds, rems);
        history.push_back(master);
    }
    
    std::pair<double, double> memory_usage() const {
        return {mvq::global_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0), 0.0}; 
    }
};
