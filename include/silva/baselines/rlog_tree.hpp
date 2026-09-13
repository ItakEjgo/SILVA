#pragma once
#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <atomic>
#include <queue>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <silva/geo/point.hpp>
#include <silva/geo/operations.hpp>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;

inline std::atomic<size_t> boost_live_mem(0);

template <typename T>
class TrackingAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    TrackingAllocator() = default;
    template <typename U> TrackingAllocator(const TrackingAllocator<U>&) {}
    template<class U> struct rebind { typedef TrackingAllocator<U> other; };
    
    T* allocate(std::size_t n) {
        boost_live_mem.fetch_add(n * sizeof(T), std::memory_order_relaxed);
        return static_cast<T*>(::operator new(n * sizeof(T)));
    }
    void deallocate(T* p, std::size_t n) {
        boost_live_mem.fetch_sub(n * sizeof(T), std::memory_order_relaxed);
        ::operator delete(p);
    }
};

template <typename T, typename U>
bool operator==(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return true; }
template <typename T, typename U>
bool operator!=(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return false; }

typedef bgi::rtree<Value, bgi::quadratic<32>, bgi::indexable<Value>, bgi::equal_to<Value>, TrackingAllocator<Value>> RTree;

namespace Rlog {

struct VersionNode {
    std::shared_ptr<VersionNode> parent;
    size_t depth;
    std::vector<Value> delta_inserts;
    std::vector<Value> delta_removes;
    std::shared_ptr<RTree> cached_base;
    size_t cumulative_log_size;
    size_t current_base_size;
    
    VersionNode() : depth(0), cumulative_log_size(0), current_base_size(0) {}
    
    ~VersionNode() {
        size_t bytes = (delta_inserts.capacity() + delta_removes.capacity()) * sizeof(Value);
        boost_live_mem.fetch_sub(bytes, std::memory_order_relaxed);
    }
    
    void add_memory(size_t bytes) {
        boost_live_mem.fetch_add(bytes, std::memory_order_relaxed);
    }
};

inline void get_query_view(std::shared_ptr<VersionNode> node, 
                           std::shared_ptr<RTree>& out_base, 
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

inline std::shared_ptr<VersionNode> map_init(const std::vector<Value>& base_points) {
    auto node = std::make_shared<VersionNode>();
    node->cached_base = std::make_shared<RTree>(base_points.begin(), base_points.end());
    node->current_base_size = base_points.size();
    return node;
}

inline void check_compact(std::shared_ptr<VersionNode> node, double p) {
    if (node->cumulative_log_size >= p * node->current_base_size) {
        std::shared_ptr<RTree> base;
        std::unordered_set<size_t> removed_ids;
        std::vector<Value> pending_inserts;
        get_query_view(node, base, removed_ids, pending_inserts);
        
        std::vector<Value> next_pts;
        next_pts.reserve(base->size() + pending_inserts.size());
        
        auto is_alive = [&](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
        std::copy_if(base->begin(), base->end(), std::back_inserter(next_pts), is_alive);
        next_pts.insert(next_pts.end(), pending_inserts.begin(), pending_inserts.end());
        
        node->cached_base = std::make_shared<RTree>(next_pts.begin(), next_pts.end());
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

struct MaxHeapCmp {
    bool operator()(const std::pair<double, Value>& a, const std::pair<double, Value>& b) const {
        if (a.first != b.first) return a.first < b.first;
        return a.second.second < b.second.second;
    }
};

inline std::vector<Value> range_report(std::shared_ptr<VersionNode> node, const geobase::Bounding_Box& q) {
    std::shared_ptr<RTree> base;
    std::unordered_set<size_t> removed_ids;
    std::vector<Value> pending_inserts;
    get_query_view(node, base, removed_ids, pending_inserts);
    
    std::vector<Value> result;
    bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
    auto is_alive = [&](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
    
    if (base) {
        base->query(bgi::intersects(box) && bgi::satisfies(is_alive), std::back_inserter(result));
    }
    
    for (const auto& val : pending_inserts) {
        if (bg::intersects(val.first, box)) {
            result.push_back(val);
        }
    }
    return result;
}

inline std::vector<Value> knn_report(std::shared_ptr<VersionNode> node, const geobase::Point& q, size_t k) {
    std::shared_ptr<RTree> base;
    std::unordered_set<size_t> removed_ids;
    std::vector<Value> pending_inserts;
    get_query_view(node, base, removed_ids, pending_inserts);
    
    auto calc_sqr_dist = [](const geobase::Point& p1, const BoostPoint& p2) {
        double dx = p1.x - p2.get<0>(), dy = p1.y - p2.get<1>();
        return dx*dx + dy*dy;
    };
    
    std::priority_queue<std::pair<double, Value>, std::vector<std::pair<double, Value>>, MaxHeapCmp> max_heap;
    for (const auto& val : pending_inserts) {
        double dist = calc_sqr_dist(q, val.first);
        if (max_heap.size() < k) max_heap.push({dist, val});
        else if (dist < max_heap.top().first || (dist == max_heap.top().first && val.second > max_heap.top().second.second)) {
            max_heap.pop(); max_heap.push({dist, val});
        }
    }
    
    if (base) {
        BoostPoint bg_q(q.x, q.y);
        auto is_alive = [&](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
        for (auto it = base->qbegin(bgi::nearest(bg_q, (unsigned)k) && bgi::satisfies(is_alive)); it != base->qend(); ++it) {
            double dist = calc_sqr_dist(q, it->first);
            if (max_heap.size() < k) {
                max_heap.push({dist, *it});
            } else if (dist < max_heap.top().first || (dist == max_heap.top().first && it->second > max_heap.top().second.second)) { 
                max_heap.pop(); max_heap.push({dist, *it}); 
            } else {
                break;
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

} // namespace Rlog