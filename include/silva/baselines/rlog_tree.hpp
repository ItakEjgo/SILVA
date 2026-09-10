#pragma once
#include <vector>
#include <memory>
#include <unordered_set>
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

inline size_t unordered_set_mem_estimate(size_t bucket_count) {
    return bucket_count * sizeof(void*); // simplistic estimate for unordered_set buckets
}

template <typename T>
class TrackingAllocator {
public:
    typedef T value_type;
    TrackingAllocator() = default;
    template <typename U> TrackingAllocator(const TrackingAllocator<U>&) {}
    
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
struct RlogBranch {
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
    shared_ptr<RTree> base_snapshot;
};

class RlogTree {
public:
    struct MaxHeapCmp {
        bool operator()(const std::pair<double, Value>& a, const std::pair<double, Value>& b) const {
            return a.first < b.first;
        }
    };
    int compaction_years;
    int last_compact_year;
    shared_ptr<RTree> snapshot;
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
    mutable size_t query_lookup_count = 0;
    mutable std::unordered_set<size_t> removed_ids;
    mutable bool cache_valid = false;
    
    void update_cache() const {
        if (!cache_valid) {
            removed_ids.clear();
            removed_ids.reserve(remove_log.size()); // 优化1：防 rehash 扩容
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            cache_valid = true;
        }
    }
    
    RlogTree() : compaction_years(99), last_compact_year(0) { snapshot = make_shared<RTree>(); }
    RlogTree(int c_years) : compaction_years(c_years), last_compact_year(0) { snapshot = make_shared<RTree>(); }

    void build_base(const std::vector<Value>& base_data) {
        snapshot = make_shared<RTree>(base_data.begin(), base_data.end());
    }
    void compact() {
        if (insert_log.empty() && remove_log.empty()) return;
        
        std::vector<Value> next_pts;
        next_pts.reserve(snapshot->size() + insert_log.size());
        
        if (!remove_log.empty()) {
            update_cache();
            auto is_alive = [this](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
            
            // 直接在从树里取点的迭代器层级，使用 is_alive 过滤 (std::copy_if)
            std::copy_if(snapshot->begin(), snapshot->end(), std::back_inserter(next_pts), is_alive);
            
            // 同样对 insert_log 进行过滤提取
            std::copy_if(insert_log.begin(), insert_log.end(), std::back_inserter(next_pts), is_alive);
        } else {
            next_pts.assign(snapshot->begin(), snapshot->end());
            next_pts.insert(next_pts.end(), insert_log.begin(), insert_log.end());
        }
        
        snapshot = make_shared<RTree>(next_pts.begin(), next_pts.end());
        
        // 优化2：真正释放内存，防止假性 OOM 和内存指标虚高
        std::vector<Value>().swap(insert_log); 
        std::vector<Value>().swap(remove_log);
        removed_ids = std::unordered_set<size_t>(); // 释放哈希表底层的 bucket 内存
        cache_valid = true; // 空表即为有效状态
    }
    void commit_inserts(const std::vector<Value>& new_pts) {
        insert_log.insert(insert_log.end(), new_pts.begin(), new_pts.end());
        
    }
    void merge(const RlogBranch& branch) {
        // 优化3：只插入数据时，根本不需要废弃关于“死亡名单”的哈希表缓存！
        if (!branch.insert_log.empty()) { insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end()); }
        if (!branch.remove_log.empty()) { remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end()); cache_valid = false; }
    }
    void check_and_compact(int current_year) {
        if (last_compact_year == 0) last_compact_year = current_year;
        if (current_year - last_compact_year >= compaction_years) {
            compact();
            last_compact_year = current_year;
        }
    }
    std::vector<Value> range_report(const geobase::Bounding_Box& q) const {
        query_lookup_count = 0;
        update_cache();
        std::vector<Value> result;
        bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
        auto is_alive = [this](Value const& v) { query_lookup_count++; return removed_ids.find(v.second) == removed_ids.end(); };
        
        // 优雅：直接在树的遍历底层完成过滤，连中间临时数组 snap_res 都省了
        snapshot->query(bgi::intersects(box) && bgi::satisfies(is_alive), std::back_inserter(result));
        
        for (const auto& val : insert_log) {
            query_lookup_count++;
            if (removed_ids.find(val.second) == removed_ids.end()) {
                // 优化4：修正语义不一致。RTree 用的是 intersects(包含边界), 这里之前用 within(不包含边界) 是有 Bug 的。
                if (bg::intersects(val.first, box)) result.push_back(val);
            }
        }
        return result;
    }
    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
        query_lookup_count = 0;
        update_cache();
        auto calc_sqr_dist = [](const geobase::Point& p1, const BoostPoint& p2) {
            double dx = p1.x - p2.get<0>(), dy = p1.y - p2.get<1>();
            return dx*dx + dy*dy;
        };
        std::priority_queue<std::pair<double, Value>, std::vector<std::pair<double, Value>>, MaxHeapCmp> max_heap;
        for (const auto& val : insert_log) {
            query_lookup_count++;
            if (removed_ids.find(val.second) == removed_ids.end()) {
                double dist = calc_sqr_dist(q, val.first);
                if (max_heap.size() < k) max_heap.push({dist, val});
                else if (dist < max_heap.top().first) { max_heap.pop(); max_heap.push({dist, val}); }
            }
        }
        
        BoostPoint bg_q(q.x, q.y);
        auto is_alive = [this](Value const& v) { query_lookup_count++; return removed_ids.find(v.second) == removed_ids.end(); };
        
        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned)k) && bgi::satisfies(is_alive)); it != snapshot->qend(); ++it) {
            double dist = calc_sqr_dist(q, it->first);
            if (max_heap.size() < k) {
                max_heap.push({dist, *it});
            } else if (dist < max_heap.top().first) { 
                max_heap.pop(); max_heap.push({dist, *it}); 
            } else {
                break; // 优化5：提前终止！如果树里找出的点已经比 max_heap 里的点更远，后续的树节点只会更远，直接打断 RTree 遍历！
            }
        }
        
        std::vector<Value> result;
        while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
        return result;
    }
    size_t size() const { return snapshot->size() + insert_log.size() + remove_log.size(); }
};
