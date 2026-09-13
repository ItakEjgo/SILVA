#pragma once
#include <string>
#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <cpam/cpam.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <cpdd/cpdd.h>

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef std::pair<BoostPoint, size_t> Value;

namespace PKDLog {

using point_t = cpdd::PointID<double, 2>;
using tree_t = cpdd::ParallelKDtree<point_t>;

inline size_t calculate_tree_memory(tree_t::node* T) {
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

struct VersionNode {
    std::shared_ptr<VersionNode> parent;
    size_t depth;
    std::vector<Value> delta_inserts;
    std::vector<Value> delta_removes;
    std::shared_ptr<tree_t> cached_base;
    size_t cumulative_log_size;
    size_t current_base_size;
    
    VersionNode() : depth(0), cumulative_log_size(0), current_base_size(0) {}
    
    ~VersionNode() {
        size_t bytes = (delta_inserts.capacity() + delta_removes.capacity()) * sizeof(Value);
        cpam::cpam_live_mem.fetch_sub(bytes, std::memory_order_relaxed);
    }
    
    void add_memory(size_t bytes) {
        cpam::cpam_live_mem.fetch_add(bytes, std::memory_order_relaxed);
    }
};

inline void get_query_view(std::shared_ptr<VersionNode> node, 
                           std::shared_ptr<tree_t>& out_base, 
                           std::unordered_set<size_t>& out_removed_ids, 
                           std::vector<Value>& out_pending_inserts) {
    std::vector<std::shared_ptr<VersionNode>> path;
    auto curr = node;
    while (curr) {
        if (curr->cached_base) {
            out_base = curr->cached_base;
            break;
        }
        path.push_back(curr);
        curr = curr->parent;
    }
    
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        for (const auto& rm : (*it)->delta_removes) {
            out_removed_ids.insert(rm.second);
        }
    }
    
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        for (const auto& ins : (*it)->delta_inserts) {
            if (out_removed_ids.find(ins.second) == out_removed_ids.end()) {
                out_pending_inserts.push_back(ins);
            }
        }
    }
}

inline std::shared_ptr<VersionNode> map_init(const std::vector<Value>& P_base) {
    auto node = std::make_shared<VersionNode>();
    
    parlay::sequence<point_t> P(P_base.size());
    parlay::parallel_for(0, P_base.size(), [&](size_t i) {
        std::array<double, 2> coords = {P_base[i].first.get<0>(), P_base[i].first.get<1>()};
        P[i] = point_t(coords, P_base[i].second);
    });
    
    auto raw_tree = new tree_t();
    raw_tree->build(parlay::make_slice(P), 2);
    size_t tree_mem = calculate_tree_memory(raw_tree->get_root());
    cpam::cpam_live_mem.fetch_add(tree_mem, std::memory_order_relaxed);
    std::shared_ptr<tree_t> tree(raw_tree, [tree_mem](tree_t* p) {
        cpam::cpam_live_mem.fetch_sub(tree_mem, std::memory_order_relaxed);
        p->delete_tree();
        delete p;
    });
    
    node->cached_base = tree;
    node->current_base_size = P_base.size();
    return node;
}

inline void check_compact(std::shared_ptr<VersionNode> node, double p) {
    if (node->cumulative_log_size >= p * node->current_base_size) {
        std::shared_ptr<tree_t> base;
        std::unordered_set<size_t> removed_ids;
        std::vector<Value> pending_inserts;
        get_query_view(node, base, removed_ids, pending_inserts);
        
        size_t tree_sz = (base && base->get_root() != nullptr) ? base->get_root()->size : 0;
        
        parlay::sequence<point_t> tree_pts;
        if (tree_sz > 0) {
            tree_pts = parlay::sequence<point_t>::uninitialized(tree_sz);
            base->flatten(base->get_root(), parlay::make_slice(tree_pts));
        }

        auto is_alive = [&](const point_t& pt) { 
            return removed_ids.find(pt.id) == removed_ids.end(); 
        };
        auto live_tree_pts = parlay::filter(tree_pts, is_alive);
        
        parlay::sequence<point_t> log_pts(pending_inserts.size());
        parlay::parallel_for(0, pending_inserts.size(), [&](size_t i) {
            std::array<double, 2> coords = {pending_inserts[i].first.get<0>(), pending_inserts[i].first.get<1>()};
            log_pts[i] = point_t(coords, pending_inserts[i].second);
        });

        parlay::sequence<point_t> next_pts = parlay::sequence<point_t>::uninitialized(live_tree_pts.size() + log_pts.size());
        parlay::parallel_for(0, live_tree_pts.size(), [&](size_t i) {
            next_pts[i] = live_tree_pts[i];
        });
        parlay::parallel_for(0, log_pts.size(), [&](size_t i) {
            next_pts[live_tree_pts.size() + i] = log_pts[i];
        });

        auto raw_tree = new tree_t();
        raw_tree->build(parlay::make_slice(next_pts), 2);
        size_t tree_mem = calculate_tree_memory(raw_tree->get_root());
        cpam::cpam_live_mem.fetch_add(tree_mem, std::memory_order_relaxed);
        std::shared_ptr<tree_t> new_tree(raw_tree, [tree_mem](tree_t* p) {
            cpam::cpam_live_mem.fetch_sub(tree_mem, std::memory_order_relaxed);
            p->delete_tree();
            delete p;
        });
        
        node->cached_base = new_tree;
        node->current_base_size = next_pts.size();
        node->cumulative_log_size = 0;
    }
}

inline std::shared_ptr<VersionNode> map_insert(std::shared_ptr<VersionNode> parent, const std::vector<Value>& insert_pts, double p = 1.0) {
    auto node = std::make_shared<VersionNode>();
    node->parent = parent;
    if (parent) {
        node->depth = parent->depth + 1;
        node->current_base_size = parent->current_base_size;
        node->cumulative_log_size = parent->cumulative_log_size + insert_pts.size();
    }
    node->delta_inserts = insert_pts;
    node->add_memory(node->delta_inserts.capacity() * sizeof(Value));
    check_compact(node, p);
    return node;
}

inline std::shared_ptr<VersionNode> map_delete(std::shared_ptr<VersionNode> parent, const std::vector<Value>& delete_pts, double p = 1.0) {
    auto node = std::make_shared<VersionNode>();
    node->parent = parent;
    if (parent) {
        node->depth = parent->depth + 1;
        node->current_base_size = parent->current_base_size;
        node->cumulative_log_size = parent->cumulative_log_size + delete_pts.size();
    }
    node->delta_removes = delete_pts;
    node->add_memory(node->delta_removes.capacity() * sizeof(Value));
    check_compact(node, p);
    return node;
}

inline std::vector<Value> range_report(std::shared_ptr<VersionNode> node, const geobase::Bounding_Box& q_copy) {
    std::shared_ptr<tree_t> base;
    std::unordered_set<size_t> removed_ids;
    std::vector<Value> pending_inserts;
    get_query_view(node, base, removed_ids, pending_inserts);
    
    std::vector<Value> result;
    auto is_alive = [&](const point_t& pt) {
        return removed_ids.find(pt.id) == removed_ids.end();
    };

    if (base && base->get_root() != nullptr) {
        tree_t::box queryBox;
        queryBox.first.pnt[0] = q_copy.first.x;
        queryBox.first.pnt[1] = q_copy.first.y;
        queryBox.second.pnt[0] = q_copy.second.x;
        queryBox.second.pnt[1] = q_copy.second.y;
        
        size_t tree_sz = base->get_root()->size;
        
        parlay::sequence<point_t> temp_out = parlay::sequence<point_t>::uninitialized(tree_sz);
        
        size_t tree_cnt = base->range_query_serial(queryBox, parlay::make_slice(temp_out));
        
        for (size_t i = 0; i < tree_cnt; i++) {
            auto& p = temp_out[i];
            if (is_alive(p)) {
                result.push_back({BoostPoint(p.pnt[0], p.pnt[1]), p.id});
            }
        }
    }
    
    bg::model::box<BoostPoint> box(BoostPoint(q_copy.first.x, q_copy.first.y), BoostPoint(q_copy.second.x, q_copy.second.y));
    for (const auto& val : pending_inserts) {
        if (bg::intersects(val.first, box)) {
            result.push_back(val);
        }
    }
    return result;
}

inline std::vector<Value> knn_report(std::shared_ptr<VersionNode> node, const geobase::Point& q, size_t k) {
    std::shared_ptr<tree_t> base;
    std::unordered_set<size_t> removed_ids;
    std::vector<Value> pending_inserts;
    get_query_view(node, base, removed_ids, pending_inserts);
    
    struct MaxHeapCmp {
        bool operator()(const std::pair<double, Value>& a, const std::pair<double, Value>& b) const {
            if (a.first != b.first) return a.first < b.first;
            return a.second.second < b.second.second;
        }
    };
    std::priority_queue<std::pair<double, Value>, std::vector<std::pair<double, Value>>, MaxHeapCmp> max_heap;
    
    auto calc_sqr_dist = [](const geobase::Point& p1, const BoostPoint& p2) {
        double dx = p1.x - p2.get<0>(), dy = p1.y - p2.get<1>();
        return dx*dx + dy*dy;
    };
    
    for (const auto& val : pending_inserts) {
        double dist = calc_sqr_dist(q, val.first);
        if (max_heap.size() < k) max_heap.push({dist, val});
        else if (dist < max_heap.top().first || (dist == max_heap.top().first && val.second > max_heap.top().second.second)) {
            max_heap.pop(); max_heap.push({dist, val});
        }
    }
    
    if (base && base->get_root() != nullptr) {
        point_t q_pt(std::array<double, 2>{q.x, q.y}, 0);
        tree_t::box nodeBox = base->get_root_box();
        using nn_pair = std::pair<point_t, double>;
        std::vector<nn_pair> out_buffer(k);
        parlay::slice<nn_pair*, nn_pair*> out_slice(out_buffer.data(), out_buffer.data() + k);
        cpdd::kBoundedQueue<point_t, nn_pair> bq(out_slice);
        
        size_t visNodeNum = 0;
        auto is_alive = [&](const point_t& pt) {
            return removed_ids.find(pt.id) == removed_ids.end();
        };
        base->k_nearest(base->get_root(), q_pt, k, bq, nodeBox, visNodeNum, is_alive);
        for (size_t i = 0; i < bq.m_count; i++) {
            auto item = out_buffer[i];
            Value v = {BoostPoint(item.first.pnt[0], item.first.pnt[1]), item.first.id};
            double dist = item.second;
            if (max_heap.size() < k) max_heap.push({dist, v});
            else if (dist < max_heap.top().first || (dist == max_heap.top().first && v.second > max_heap.top().second.second)) {
                max_heap.pop(); max_heap.push({dist, v});
            }
        }
    }
    
    std::vector<Value> result;
    while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
    return result;
}

struct diff_type {
    std::vector<Value> add;
    std::vector<Value> remove;
    void compact() {}
};

inline void map_spatial_diff(std::shared_ptr<VersionNode> v1, std::shared_ptr<VersionNode> v2, const geobase::Bounding_Box& q, diff_type& out_diff) {
    if (v1 == v2) return;
    
    std::shared_ptr<VersionNode> ptr1 = v1;
    std::shared_ptr<VersionNode> ptr2 = v2;
    
    while (ptr1 && ptr2 && ptr1->depth > ptr2->depth) ptr1 = ptr1->parent;
    while (ptr1 && ptr2 && ptr2->depth > ptr1->depth) ptr2 = ptr2->parent;
    
    while (ptr1 && ptr2 && ptr1 != ptr2) {
        ptr1 = ptr1->parent;
        ptr2 = ptr2->parent;
    }
    auto lca = ptr1;
    
    std::unordered_map<size_t, int> state1;
    std::unordered_map<size_t, int> state2;
    std::unordered_map<size_t, Value> values;
    
    auto trace_branch = [&](std::shared_ptr<VersionNode> leaf, std::unordered_map<size_t, int>& state) {
        auto curr = leaf;
        while (curr && curr != lca) {
            for (const auto& rm : curr->delta_removes) {
                if (state.find(rm.second) == state.end()) state[rm.second] = -1;
                values[rm.second] = rm;
            }
            for (const auto& ins : curr->delta_inserts) {
                if (state.find(ins.second) == state.end()) state[ins.second] = 1;
                else if (state[ins.second] == -1) state[ins.second] = 0;
                values[ins.second] = ins;
            }
            curr = curr->parent;
        }
    };
    
    trace_branch(v1, state1);
    trace_branch(v2, state2);
    
    bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
    
    for (const auto& kv : values) {
        size_t id = kv.first;
        if (!bg::intersects(kv.second.first, box)) continue;
        
        int s1 = state1.count(id) ? state1[id] : 0;
        int s2 = state2.count(id) ? state2[id] : 0;
        
        int net_change = s2 - s1;
        if (net_change > 0) {
            out_diff.add.push_back(kv.second);
        } else if (net_change < 0) {
            out_diff.remove.push_back(kv.second);
        }
    }
}

} // namespace PKDLog