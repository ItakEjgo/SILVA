#pragma once
#include <string>
#include <vector>
#include "runner_common.hpp"

// Using the same tracking allocator and types from RlogTree
typedef bgi::rtree<Value, bgi::quadratic<32>, bgi::indexable<Value>, bgi::equal_to<Value>, TrackingAllocator<Value>> RawBoostRTree;

struct BoostRunner {
    static std::string name() { return "BoostRtree"; }
    
    RawBoostRTree* tree;
    
    BoostRunner() { tree = new RawBoostRTree(); }
    ~BoostRunner() { delete tree; }
    
    void build_base(const std::vector<Value>& P_base) {
        delete tree;
        tree = new RawBoostRTree(P_base.begin(), P_base.end());
    }
    
    void reset_query_stats() {}
    
    size_t get_nodes_touched() const { return 0; }
    
    void range_query(geobase::Bounding_Box q_copy, parlay::sequence<geobase::Point>& /*shared_out*/, size_t& cnt, size_t& h) {
        bg::model::box<BoostPoint> box(BoostPoint(q_copy.first.x, q_copy.first.y), BoostPoint(q_copy.second.x, q_copy.second.y));
        std::vector<Value> result;
        tree->query(bgi::intersects(box), std::back_inserter(result));
        cnt = result.size();
        for(const auto& v : result) h += v.second;
    }
    
    void knn_query(geobase::Point q, size_t k, size_t& cnt, size_t& h) {
        BoostPoint bg_q(q.x, q.y);
        std::vector<Value> result;
        tree->query(bgi::nearest(bg_q, (unsigned)k), std::back_inserter(result));
        cnt = result.size();
        for(const auto& v : result) h += v.second;
    }
    
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        for (const auto& op : rems) {
            tree->remove(std::make_pair(BoostPoint(op.x, op.y), op.id));
        }
        for (const auto& op : adds) {
            tree->insert(std::make_pair(BoostPoint(op.x, op.y), op.id));
        }
    }
    
    std::pair<double, double> memory_usage() const {
        double rtree_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
        return {rtree_mem, 0.0};
    }
};
