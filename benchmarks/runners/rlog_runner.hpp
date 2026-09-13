#pragma once
#include <string>
#include <vector>
#include "runner_common.hpp"
#include <silva/baselines/rlog_tree.hpp>

struct RlogRunner {
    static std::string name() { return "RLogTree"; }
    
    std::shared_ptr<Rlog::VersionNode> tree;
    
    RlogRunner(int /*compaction_years*/ = 1) {}
    
    ~RlogRunner() {}
    
    void build_base(const std::vector<Value>& P_base) {
        tree = Rlog::map_init(P_base);
    }
    
    void reset_query_stats() {
        // Not implemented in functional version
    }
    
    size_t get_nodes_touched() const {
        return 0; // Not implemented
    }
    
    void range_query(const geobase::Bounding_Box& q_copy, parlay::sequence<geobase::Point>& /*shared_out*/, size_t& cnt, size_t& h) {
        auto res = Rlog::range_report(tree, q_copy);
        cnt = res.size();
        for(const auto& v : res) h += v.second;
    }
    
    void knn_query(geobase::Point q, size_t k, size_t& cnt, size_t& h) {
        auto res = Rlog::knn_report(tree, q, k);
        cnt = res.size();
        for(const auto& v : res) h += v.second;
    }
    
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        std::vector<Value> boost_adds;
        for(const auto& op : adds) boost_adds.push_back(std::make_pair(BoostPoint(op.x, op.y), op.id));
        
        std::vector<Value> boost_rems;
        for(const auto& r_pt : rems) boost_rems.push_back(std::make_pair(BoostPoint(r_pt.x, r_pt.y), r_pt.id));
        
        tree = Rlog::map_insert(tree, boost_adds, 0.2); // p=0.2
        tree = Rlog::map_delete(tree, boost_rems, 0.2);
    }
    
    std::pair<double, double> memory_usage() const {
        double mem_mb = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
        return {mem_mb, mem_mb};
    }
};