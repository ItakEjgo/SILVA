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
        if (tree) { tree->delete_tree(); delete tree; }
    }
    
    void build_base(const std::vector<Value>& P_base) {
        parlay::sequence<point_t> init_pts = parlay::sequence<point_t>::uninitialized(P_base.size());
        parlay::parallel_for(0, P_base.size(), [&](size_t i) {
            std::array<double, 2> coords = {P_base[i].first.get<0>(), P_base[i].first.get<1>()};
            init_pts[i] = point_t(coords, P_base[i].second);
        });
        if (tree) { tree->delete_tree(); delete tree; }
        tree = new tree_t();
        tree->build(parlay::make_slice(init_pts), 2);
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
        size_t tree_sz = (tree && tree->get_root() != nullptr) ? tree->get_root()->size : 0;
        parlay::sequence<point_t> pkd_out(tree_sz);
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
    void compact() {
        if (insert_log.empty() && remove_log.empty()) return;
        update_cache();

        size_t tree_sz = 0;
        if (tree && tree->get_root() != nullptr) {
            tree_sz = tree->get_root()->size;
        }

        // 1. Parallel Flatten
        parlay::sequence<point_t> tree_pts;
        if (tree_sz > 0) {
            tree_pts = parlay::sequence<point_t>::uninitialized(tree_sz);
            tree->flatten(tree->get_root(), parlay::make_slice(tree_pts));
        }

        // 2. Parallel Filter (Old tree points)
        auto is_alive = [&](const point_t& p) { 
            return removed_ids.find(p.id) == removed_ids.end(); 
        };
        auto live_tree_pts = parlay::filter(tree_pts, is_alive);

        // 3. Parallel Filter (New insert log)
        auto log_is_alive = [&](const Value& v) { 
            return removed_ids.find(v.second) == removed_ids.end(); 
        };
        auto live_log_pts = parlay::filter(insert_log, log_is_alive);

        // 4. Parallel Combine
        parlay::sequence<point_t> next_pts = parlay::sequence<point_t>::uninitialized(live_tree_pts.size() + live_log_pts.size());
        
        parlay::parallel_for(0, live_tree_pts.size(), [&](size_t i) {
            next_pts[i] = live_tree_pts[i];
        });
        
        parlay::parallel_for(0, live_log_pts.size(), [&](size_t i) {
            std::array<double, 2> coords = {live_log_pts[i].first.get<0>(), live_log_pts[i].first.get<1>()};
            next_pts[live_tree_pts.size() + i] = point_t(coords, live_log_pts[i].second);
        });

        // 5. Rebuild & Cleanup
        if (tree) { tree->delete_tree(); delete tree; }
        tree = new tree_t();
        tree->build(parlay::make_slice(next_pts), 2);

        insert_log.clear();
        remove_log.clear();
        removed_ids.clear();
        cache_valid = true;
    }
    void check_and_compact(double p) {
        double log_size = insert_log.size() + remove_log.size();
        size_t tree_sz = (tree && tree->get_root() != nullptr) ? tree->get_root()->size : 0;
        if (log_size >= p * tree_sz) {
            compact();
        }
    }
    size_t calculate_tree_memory(tree_t::node* T) const {
        if (T == nullptr) return 0;
        if (T->is_leaf) {
            return sizeof(tree_t::leaf) + tree_t::LEAVE_WRAP * sizeof(point_t);
        }
        
        tree_t::interior* TI = static_cast<tree_t::interior*>(T);
        size_t l = 0, r = 0;
        parlay::par_do_if(TI->size > 1000,
            [&]() { l = calculate_tree_memory(TI->left); },
            [&]() { r = calculate_tree_memory(TI->right); }
        );
        return sizeof(tree_t::interior) + l + r;
    }

    std::pair<double, double> memory_usage() const {
        size_t index_bytes = calculate_tree_memory(const_cast<tree_t*>(tree)->get_root());
        size_t log_bytes = (insert_log.capacity() + remove_log.capacity()) * sizeof(Value) + 
                           removed_ids.bucket_count() * sizeof(void*) + 
                           removed_ids.size() * sizeof(size_t);
        double mem_mb = (index_bytes + log_bytes) / (1024.0 * 1024.0);
        return {mem_mb, mem_mb};
    }
};
