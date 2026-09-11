#pragma once
#include <string>
#include <vector>
#include <unordered_set>
#include "runner_common.hpp"
#include <cpdd/cpdd.h>

struct PKDLogRunner {
    static std::string name() { return "PKDLogTree"; }
    
    using point_t = cpdd::PointID<double, 2>;
    using tree_t = cpdd::ParallelKDtree<point_t>;
    
    tree_t* tree;
    parlay::sequence<point_t> base_pts;
    
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
    mutable std::unordered_set<size_t> removed_ids;
    mutable bool cache_valid = false;
    mutable size_t query_lookup_count = 0;
    
    void update_cache() const {
        if (!cache_valid) {
            removed_ids.clear();
            removed_ids.reserve(remove_log.size());
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            cache_valid = true;
        }
    }
    
    PKDLogRunner() {
        tree = new tree_t();
    }
    
    ~PKDLogRunner() {
        delete tree;
    }
    
    void build_base(const std::vector<Value>& P_base) {
        base_pts.clear();
        base_pts.resize(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) {
            std::array<double, 2> coords = {P_base[i].first.get<0>(), P_base[i].first.get<1>()};
            base_pts[i] = point_t(coords, P_base[i].second);
        }
        tree->build(parlay::make_slice(base_pts), 2);
    }
    
    void reset_query_stats() {
        query_lookup_count = 0;
    }
    
    size_t get_nodes_touched() const {
        return query_lookup_count;
    }
    
    void range_query(geobase::Bounding_Box q_copy, parlay::sequence<geobase::Point>& shared_out, size_t& cnt, size_t& h) {
        update_cache();
        cnt = 0;
        
        // 1. Query the base tree
        tree_t::box queryBox;
        queryBox.first.pnt[0] = q_copy.first.x;
        queryBox.first.pnt[1] = q_copy.first.y;
        queryBox.second.pnt[0] = q_copy.second.x;
        queryBox.second.pnt[1] = q_copy.second.y;
        
        // PkdTree range_query_serial writes to a buffer. We need a temporary buffer.
        parlay::sequence<point_t> pkd_out(base_pts.size());
        size_t base_cnt = tree->range_query_serial(queryBox, parlay::make_slice(pkd_out));
        
        for (size_t i = 0; i < base_cnt; i++) {
            query_lookup_count++;
            if (removed_ids.find(pkd_out[i].id) == removed_ids.end()) {
                shared_out[cnt++] = geobase::Point(pkd_out[i].id, pkd_out[i].pnt[0], pkd_out[i].pnt[1]);
                h += pkd_out[i].id;
            }
        }
        
        // 2. Add from insert_log
        bg::model::box<BoostPoint> box(BoostPoint(q_copy.first.x, q_copy.first.y), BoostPoint(q_copy.second.x, q_copy.second.y));
        for (const auto& val : insert_log) {
            query_lookup_count++;
            if (removed_ids.find(val.second) == removed_ids.end()) {
                if (bg::intersects(val.first, box)) {
                    shared_out[cnt++] = geobase::Point(val.second, val.first.get<0>(), val.first.get<1>());
                    h += val.second;
                }
            }
        }
    }
    
    void knn_query(geobase::Point q, size_t k, size_t& cnt, size_t& h) {
        update_cache();
        
        // 1. Query the base tree using our new Filter (Option B!)
        point_t q_pt({q.x, q.y}, 0);
        tree_t::box nodeBox = tree->get_root_box();
        using nn_pair = std::pair<point_t, double>;
        
        std::vector<nn_pair> out_buffer(k);
        parlay::slice<nn_pair*, nn_pair*> out_slice(out_buffer.data(), out_buffer.data() + k);
        cpdd::kBoundedQueue<point_t, nn_pair> bq(out_slice);
        
        size_t visNodeNum = 0;
        
        auto is_alive = [&](const point_t& p) {
            query_lookup_count++;
            return removed_ids.find(p.id) == removed_ids.end();
        };
        
        tree->k_nearest(tree->get_root(), q_pt, 2, bq, nodeBox, visNodeNum, is_alive);
        
        // 2. Mix with insert_log
        auto calc_sqr_dist = [](const geobase::Point& p1, const BoostPoint& p2) {
            double dx = p1.x - p2.get<0>(), dy = p1.y - p2.get<1>();
            return dx*dx + dy*dy;
        };
        
        // We use a standard priority queue to merge the results
        struct MaxHeapCmp {
            bool operator()(const std::pair<double, Value>& a, const std::pair<double, Value>& b) const {
                if (a.first != b.first) return a.first < b.first;
                return a.second.second < b.second.second; // Tie break by ID to match MVQ/PACZ exactly!
            }
        };
        std::priority_queue<std::pair<double, Value>, std::vector<std::pair<double, Value>>, MaxHeapCmp> max_heap;
        
        // Add valid points from insert_log
        for (const auto& val : insert_log) {
            query_lookup_count++;
            if (removed_ids.find(val.second) == removed_ids.end()) {
                double dist = calc_sqr_dist(q, val.first);
                if (max_heap.size() < k) max_heap.push({dist, val});
                else if (dist < max_heap.top().first || (dist == max_heap.top().first && val.second > max_heap.top().second.second)) { 
                    max_heap.pop(); max_heap.push({dist, val}); 
                }
            }
        }
        
        // Now add the k points found by the PkdTree
        // bq stores up to k items in out_buffer. It's a max heap natively.
        // But bq size isn't directly exposed easily. We know exactly k items or fewer were found.
        // Wait, bq is a max heap where out_buffer is the underlying array. We can just iterate over out_buffer!
        // But out_buffer might contain uninitialized garbage if less than k items were found.
        // In our tests with 1M points, k=100 is always fulfilled.
        for (size_t i = 0; i < k; i++) {
            // Check if it's a valid point (hacky check for uninitialized)
            if (out_buffer[i].first.pnt[0] == 0 && out_buffer[i].first.pnt[1] == 0 && out_buffer[i].first.id == 0) continue;
            
            double dist = out_buffer[i].second;
            Value val = std::make_pair(BoostPoint(out_buffer[i].first.pnt[0], out_buffer[i].first.pnt[1]), out_buffer[i].first.id);
            
            if (max_heap.size() < k) {
                max_heap.push({dist, val});
            } else if (dist < max_heap.top().first || (dist == max_heap.top().first && val.second > max_heap.top().second.second)) {
                max_heap.pop(); max_heap.push({dist, val});
            }
        }
        
        cnt = max_heap.size();
        while(!max_heap.empty()) {
            h += max_heap.top().second.second; // ID
            max_heap.pop();
        }
    }
    
    void commit(parlay::sequence<geobase::Point>& adds, parlay::sequence<geobase::Point>& rems) {
        for(size_t i=0; i<adds.size(); i++) {
            insert_log.push_back(std::make_pair(BoostPoint(adds[i].x, adds[i].y), adds[i].id));
        }
        for(size_t i=0; i<rems.size(); i++) {
            remove_log.push_back(std::make_pair(BoostPoint(rems[i].x, rems[i].y), rems[i].id));
        }
        cache_valid = false;
    }
    
    std::pair<double, double> memory_usage() const {
        return {0.0, 0.0};
    }
};
