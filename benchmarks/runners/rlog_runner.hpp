#pragma once
#include <string>
#include <vector>
#include "runner_common.hpp"
#include <silva/baselines/rlog_tree.hpp>
#include <silva/core/benchmark_utils.hpp>

struct RlogRunner {
    static std::string name() { return "RLogTree"; }
    
    RlogTree* tree;
    std::vector<RlogBranch> history;
    
    RlogRunner(int compaction_years) {
        tree = new RlogTree(compaction_years);
    }
    
    ~RlogRunner() {
        delete tree;
    }
    
    void build_base(const std::vector<Value>& P_base) {
        std::vector<Value> P_base_conv(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) {
            P_base_conv[i] = std::make_pair(BoostPoint(P_base[i].first.get<0>(), P_base[i].first.get<1>()), P_base[i].second);
        }
        tree->build_base(P_base_conv);
    }
    
    void reset_query_stats() {
        tree->query_lookup_count = 0;
    }
    
    size_t get_nodes_touched() const {
        return tree->query_lookup_count;
    }
    
    void range_query(const geobase::Bounding_Box& q_copy, parlay::sequence<geobase::Point>& /*shared_out*/, size_t& cnt, size_t& h) {
        auto res = tree->range_report(q_copy);
        cnt = res.size();
        for(const auto& v : res) h += v.second;
    }
    
    void knn_query(geobase::Point q, size_t k, size_t& cnt, size_t& h) {
        auto res = tree->knn_report(q, k);
        cnt = res.size();
        for(const auto& v : res) h += v.second;
    }
    
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        std::vector<Value> boost_adds;
        for(const auto& op : adds) boost_adds.push_back(std::make_pair(BoostPoint(op.x, op.y), op.id));
        
        RlogBranch branch;
        branch.insert_log = boost_adds;
        for(const auto& r_pt : rems) branch.remove_log.push_back(std::make_pair(BoostPoint(r_pt.x, r_pt.y), r_pt.id));
        
        tree->merge(branch);
        history.push_back(branch);
        tree->check_and_compact((int)history.size()); // Simulate year increment
    }
    
    std::pair<double, double> memory_usage() const {
        return mem_rlog(*tree, history);
    }
};
