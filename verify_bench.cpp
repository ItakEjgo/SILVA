#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <sys/stat.h>
#include <unordered_set>
#include <atomic>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "parlay/internal/get_time.h"
#include <cpam/parse_command_line.h>
#include "helper/time_loop.h"

#include "src/mvq.hpp"
#include "src/cpamz.hpp"
#include "src/global_config.hpp"
#include "src/cpambb.hpp"

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;


inline std::atomic<size_t> boost_live_mem(0);

size_t unordered_set_mem_estimate(size_t bucket_count) {
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

struct BranchRes { double fork_ms; double commit_ms; double merge_ms; double mem_mb; };








std::pair<double, double> mem_cpambb(const CPAMBB::zmap& latest_branch) {
    return {cpam::cpam_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0), 0.0};
}


std::pair<double, double> mem_mvzd() { 
    return {mvq::global_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0), 0.0}; 
}

std::pair<double, double> mem_rlog(const RlogTree& master, const std::vector<RlogBranch>& history) {
    double rtree_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
    double delta_mem = 0;
    for (const auto& b : history) {
        delta_mem += (b.insert_log.capacity() + b.remove_log.capacity()) * sizeof(Value);
    }
    delta_mem += (master.insert_log.capacity() + master.remove_log.capacity()) * sizeof(Value);
    delta_mem += unordered_set_mem_estimate(master.removed_ids.bucket_count()) + master.removed_ids.size() * sizeof(size_t);
    return {rtree_mem, delta_mem / (1024.0 * 1024.0)};
}
void log_branch(ofstream& fout, const string& algo, int batch_idx, BranchRes res) {
    fout << std::flush;
    fout << algo << " | " << batch_idx << " | " << res.fork_ms << " | " << res.commit_ms << " | " << res.merge_ms << " | " << res.mem_mb << "\n";
}

struct YearData { int year; vector<parlay::sequence<geobase::Point>> batches; };

template<typename F_Init, typename F_Run, typename F_End>
double run_batch_time_loop(F_Init initf, F_Run runf, F_End endf) {
    return time_loop(3, 1.0, initf, runf, endf) * 1000.0;
}



struct Operation {
    char op_type;
    size_t node_id;
    double lon;
    double lat;
};

struct Event {
    int type; // 0 = START, 1 = END
    string timestamp;
    size_t changeset_id;
    std::vector<Operation> ops;
    
    bool operator<(const Event& other) const {
        if (timestamp != other.timestamp) return timestamp < other.timestamp;
        return type < other.type; 
    }
};

class CommitStream {
    std::ifstream fin;
    int current_year;
    std::string line;
    bool has_next;
    size_t last_cs_id;
    std::vector<Operation> buffer_ops;

    std::string ds_dir;
    std::string commits_dir;
    int end_year;

public:
    CommitStream(int start_year, std::string dir, std::string cd, int ey) : current_year(start_year), has_next(true), last_cs_id(0), ds_dir(dir), commits_dir(cd), end_year(ey) {
        open_year(current_year);
    }
    
    void open_year(int year) {
        if (fin.is_open()) fin.close();
        fin.open(ds_dir + "/" + commits_dir + "/commits_" + to_string(year) + ".csv");
        if (fin.is_open()) {
            getline(fin, line); // skip header
        }
        // DO NOT set has_next = false here! Missing a year shouldn't terminate the whole stream.
        // The stream will naturally terminate when current_year > end_year.
    }

    bool next_changeset(Event& out_event, int& out_year) {
        out_event.ops.clear();
        if (!has_next && buffer_ops.empty()) return false;
        
        out_year = current_year;
        
        if (!buffer_ops.empty()) {
            out_event.changeset_id = last_cs_id;
            out_event.ops = std::move(buffer_ops);
            buffer_ops.clear();
        }
        
        while (true) {
            if (!getline(fin, line)) {
                current_year++;
                if (current_year > end_year) { 
                    has_next = false; 
                    return !out_event.ops.empty(); 
                }
                open_year(current_year);
                continue;
            }
            
            stringstream ss(line);
            string cs_str, start_str, end_str, op_str, nid_str, ver_str, lon_str, lat_str, ts_str;
            getline(ss, cs_str, ','); getline(ss, start_str, ','); getline(ss, end_str, ',');
            getline(ss, op_str, ','); getline(ss, nid_str, ','); getline(ss, ver_str, ',');
            getline(ss, lon_str, ','); getline(ss, lat_str, ',');
            
            size_t cs = stoull(cs_str);
            Operation op;
            op.op_type = op_str[0];
            op.node_id = stoull(nid_str);
            if (op.op_type != 'D') {
                op.lon = stod(lon_str) * 1000000.0;
                op.lat = stod(lat_str) * 1000000.0;
            } else {
                op.lon = 0; op.lat = 0;
            }
            
            if (out_event.ops.empty()) {
                out_event.changeset_id = cs;
            }
            
            if (out_event.changeset_id != cs) {
                last_cs_id = cs;
                buffer_ops.push_back(op);
                return true;
            }
            
            out_event.ops.push_back(op);
        }
    }
};

int main(int argc, char** argv) {
    srand(42);
    cpam::commandLine cmd(argc, argv, "");
    string run_algo = cmd.getOptionValue("-algo", "all");
    string dataset_dir = cmd.getOptionValue("-dir", "/data/bhuan102/SILVA-dataset/bhutan_workload");
    int just_build = cmd.getOptionIntValue("-just_build", 0);
    int start_year = cmd.getOptionIntValue("-start_year", 2018);
    int end_year = cmd.getOptionIntValue("-end_year", 2026);
    size_t q_step = cmd.getOptionIntValue("-q_step", 1000);
    string commits_dir = cmd.getOptionValue("-commits_dir", "01_commits");

    std::cout << "[Init] Pre-scanning for Dynamic MBR..." << std::endl;
    double min_x = 1e18, min_y = 1e18, max_x = -1e18, max_y = -1e18;
    
    auto update_mbr = [&](double lon, double lat) {
        if (lon < min_x) min_x = lon;
        if (lat < min_y) min_y = lat;
        if (lon > max_x) max_x = lon;
        if (lat > max_y) max_y = lat;
    };

    // Scan Base
    ifstream fin_base(dataset_dir + "/00_build/base_snapshot_" + to_string(start_year - 1) + ".csv");
    string line;
    getline(fin_base, line);
    while (getline(fin_base, line)) {
        stringstream ss(line);
        string nid_str, ver_str, lon_str, lat_str;
        getline(ss, nid_str, ','); getline(ss, ver_str, ',');
        getline(ss, lon_str, ','); getline(ss, lat_str, ',');
        update_mbr(stod(lon_str) * 1000000.0, stod(lat_str) * 1000000.0);
    }
    fin_base.close();

    // Scan Commits
    for (int year = start_year; year <= end_year; year++) {
        ifstream fin(dataset_dir + "/" + commits_dir + "/commits_" + to_string(year) + ".csv");
        if (!fin.is_open()) continue;
        getline(fin, line);
        while(getline(fin, line)) {
            stringstream ss(line);
            string cs_str, start_str, end_str, op_str, nid_str, ver_str, lon_str, lat_str;
            getline(ss, cs_str, ','); getline(ss, start_str, ','); getline(ss, end_str, ',');
            getline(ss, op_str, ','); getline(ss, nid_str, ','); getline(ss, ver_str, ',');
            getline(ss, lon_str, ','); getline(ss, lat_str, ',');
            if (op_str[0] != 'D') update_mbr(stod(lon_str) * 1000000.0, stod(lat_str) * 1000000.0);
        }
    }

    // 1% Padding
    double dx = (max_x - min_x) * 0.01;
    double dy = (max_y - min_y) * 0.01;
    mvq::Config::get().largest_mbr = geobase::Bounding_Box(
        geobase::Point(min_x - dx, min_y - dy), 
        geobase::Point(max_x + dx, max_y + dy)
    );
    std::cout << "[Init] Dynamic MBR set: X[" << min_x-dx << ", " << max_x+dx << "] Y[" << min_y-dy << ", " << max_y+dy << "]" << std::endl;

    std::cout << "[Init] Loading Base Snapshot (" << (start_year - 1) << ")..." << std::endl;
    std::vector<geobase::Point> P_base_vec;
    fin_base.open(dataset_dir + "/00_build/base_snapshot_" + to_string(start_year - 1) + ".csv");
    getline(fin_base, line);
    while (getline(fin_base, line)) {
        stringstream ss(line);
        string nid_str, ver_str, lon_str, lat_str;
        getline(ss, nid_str, ','); getline(ss, ver_str, ',');
        getline(ss, lon_str, ','); getline(ss, lat_str, ',');
        P_base_vec.push_back(geobase::Point(stoull(nid_str), stod(lon_str)*1000000.0, stod(lat_str)*1000000.0));
    }
    parlay::sequence<geobase::Point> P_base(P_base_vec.begin(), P_base_vec.end());
    auto P_base_set = geobase::get_sorted_points(P_base);

    std::unordered_map<size_t, geobase::Point> active_nodes;
    std::vector<geobase::Point> all_ever_existed_nodes;
    for(auto& pt : P_base_vec) {
        active_nodes[pt.id] = pt;
        all_ever_existed_nodes.push_back(pt);
    }

    double base_build_ms = 0.0;
    // Tree definitions
    mvq::Tree* mvzd_tree = nullptr;
    shared_ptr<mvq::BaseNode> mvzd_master;
    CPAMBB::zmap cpambb_master;
    RlogTree *rlog_master[6];
    std::vector<shared_ptr<mvq::BaseNode>> mvzd_global_history;
    std::vector<CPAMBB::zmap> cpambb_global_history;
    std::vector<RlogBranch> rlog_global_history[6];

    auto sub_pts = P_base_set.substr(0, min((size_t)10000, P_base_set.size()));
    if (run_algo == "MVZD" || run_algo == "all") {
        mvq::Tree dummy(mvq::Config::get().leaf_size);
        dummy.build(sub_pts);
        mvzd_tree = new mvq::Tree(mvq::Config::get().leaf_size);
        cpam::timer t_b;
        mvzd_tree->build(P_base_set);
        base_build_ms = t_b.stop() * 1000.0;
        mvzd_master = mvzd_tree->root;
        mvzd_global_history.push_back(mvzd_master);
    }
    if (run_algo == "CPAMBB" || run_algo == "all") {
        auto dummy = CPAMBB::map_init(sub_pts, false);
        cpam::timer t_b;
        cpambb_master = CPAMBB::map_init(P_base_set, false);
        base_build_ms = t_b.stop() * 1000.0;
        cpambb_global_history.push_back(cpambb_master);
    }
    std::vector<std::pair<BoostPoint, size_t>> P_base_conv;
    if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
        P_base_conv.resize(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) {
            P_base_conv[i] = make_pair(BoostPoint(P_base[i].x, P_base[i].y), P_base[i].id);
        }
        string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
        int thresholds[6] = { 1, 2, 3, 4, 5, 1000 };
        for(int i=0; i<6; i++) {
            if (run_algo == names[i] || run_algo == "all") {
                rlog_master[i] = new RlogTree(thresholds[i]);
                cpam::timer t_b;
                rlog_master[i]->build_base(P_base_conv);
                base_build_ms = t_b.stop() * 1000.0;
                RlogBranch init_b;
                rlog_global_history[i].push_back(init_b);
            }
        }
    }
        
    if (just_build) {
        std::cout << "BUILD_TIME_MS|" << run_algo << "|" << base_build_ms << std::endl;
        return 0;
    }
    
    
    // Core Query Lambda (Dynamic Generation embedded inside)
    auto run_hybrid_queries = [&](string fname) {
        ofstream fout(fname, std::ios_base::app);
        if (!fout.is_open()) return;
        
        size_t n = active_nodes.size();
        if (n == 0) return;
        
        // Global sort overhead completely eliminated! 
        // We now use rejection sampling on all_ever_existed_nodes which is O(1)
        srand(n); // Ensure deterministic repeatable sequences per checkpoint
        
        std::vector<size_t> K_list = {1, 10, 100};
        
        parlay::sequence<geobase::Point> current_knn(100);
        for(int i=0; i<100; i++) {
            while(true) {
                auto pt = all_ever_existed_nodes[rand() % all_ever_existed_nodes.size()];
                if (active_nodes.find(pt.id) != active_nodes.end()) {
                    current_knn[i] = active_nodes[pt.id];
                    break;
                }
            }
        }
        
        parlay::sequence<geobase::Bounding_Box> qs_small(100), qs_med(100), qs_large(100);
        
        std::cout << "    [Gen Queries] Spatial Rectangle Heuristic. N=" << n << std::endl;
        
        double world_w = mvq::Config::get().largest_mbr.second.x - mvq::Config::get().largest_mbr.first.x;
        double world_h = mvq::Config::get().largest_mbr.second.y - mvq::Config::get().largest_mbr.first.y;
        
        auto build_qs = [&](double ratio, parlay::sequence<geobase::Bounding_Box>& qs) {
            double dx = world_w * ratio;
            double dy = world_h * ratio;
            for(int i=0; i<100; i++) {
                geobase::Point center;
                while(true) {
                    auto pt = all_ever_existed_nodes[rand() % all_ever_existed_nodes.size()];
                    if (active_nodes.find(pt.id) != active_nodes.end()) {
                        center = active_nodes[pt.id];
                        break;
                    }
                }
                qs[i] = geobase::Bounding_Box(
                    geobase::Point(center.x - dx, center.y - dy),
                    geobase::Point(center.x + dx, center.y + dy)
                );
            }
        };
        
        // Small: 0.005% of width/height, Med: 0.05%, Large: 0.5% (Scaled down 10x)
        build_qs(0.00005, qs_small);
        build_qs(0.0005, qs_med);
        build_qs(0.005, qs_large);
        
        // As requested: dynamically sized std::vector strictly bounded by N to completely eliminate buffer overflows
        parlay::sequence<geobase::Point> shared_out(n + 100000);

        auto run_range_set = [&](const parlay::sequence<geobase::Bounding_Box>& qs, string q_label) {
            if (run_algo == "MVZD" || run_algo == "all") {
                for(size_t i=0; i<qs.size(); i++) {
                    geobase::Bounding_Box q_copy = qs[i];
                    size_t cnt = 0;
                    size_t h = 0;
                    cpam::timer t_q;
                    for (int rep = 0; rep < 3; rep++) {
                        cnt = 0; h = 0; mvq::query_nodes_touched = 0;
                        mvzd_tree->range_report(mvzd_master, q_copy, mvq::Config::get().largest_mbr, cnt, shared_out);
                        for(size_t j=0; j<cnt; j++) h += shared_out[j].id;
                    }
                    double q_ms = (t_q.stop() * 1000.0) / 3.0;
                    auto mems = mem_mvzd();
                    fout << "MVZD | " << q_label << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << mvq::query_nodes_touched.load() << "\n";
                }
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                for(size_t i=0; i<qs.size(); i++) {
                    geobase::Bounding_Box q_copy = qs[i];
                    size_t cnt = 0;
                    size_t h = 0;
                    cpam::timer t_q;
                    for (int rep = 0; rep < 3; rep++) {
                        cpam::cpam_query_nodes_touched = 0;
                        h = 0;
                        cnt = CPAMBB::range_report(cpambb_master, q_copy, shared_out, false);
                        for(size_t j=0; j<cnt; j++) h += shared_out[j].id;
                    }
                    double q_ms = (t_q.stop() * 1000.0) / 3.0;
                    auto mems = mem_cpambb(cpambb_master);
                    size_t est_nodes = cpam::cpam_query_nodes_touched.load();
                    fout << "CPAMBB | " << q_label << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << est_nodes << "\n";
                }
            }
            if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                for(int j=0; j<6; j++) {
                    if (run_algo == names[j] || run_algo == "all") {
                        for(size_t i=0; i<qs.size(); i++) {
                            size_t cnt = 0;
                            size_t h = 0;
                            cpam::timer t_q;
                            for (int rep = 0; rep < 3; rep++) {
                                h = 0;
                                auto res = rlog_master[j]->range_report(qs[i]);
                                cnt = res.size();
                                for(const auto& v : res) h += v.second;
                            }
                            double q_ms = (t_q.stop() * 1000.0) / 3.0;
                            size_t est_nodes = rlog_master[j]->query_lookup_count;
                            auto mems = mem_rlog(*rlog_master[j], rlog_global_history[j]);
                            fout << names[j] << " | " << q_label << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << est_nodes << "\n";
                        }
                    }
                }
            }
        };

        run_range_set(qs_small, "Range_Small");
        run_range_set(qs_med, "Range_Med");
        run_range_set(qs_large, "Range_Large");

        for(size_t K : K_list) {
            string qk = "KNN_" + to_string(K);
            if (run_algo == "MVZD" || run_algo == "all") {
                for(size_t i=0; i<100; i++) {
                    size_t cnt = 0;
                    size_t h = 0;
                    cpam::timer t_q;
                    for (int rep = 0; rep < 3; rep++) {
                        h = 0; mvq::query_nodes_touched = 0;
                        auto res = mvzd_tree->knn_report(K, current_knn[i], mvq::Config::get().largest_mbr);
                        cnt = res.size();
                        while(!res.empty()) { h += res.top().first.id; res.pop(); }
                    }
                    double q_ms = (t_q.stop() * 1000.0) / 3.0;
                    auto mems = mem_mvzd();
                    fout << "MVZD | " << qk << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << mvq::query_nodes_touched.load() << "\n";
                }
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                for(size_t i=0; i<100; i++) {
                    size_t cnt = 0;
                    size_t h = 0;
                    cpam::timer t_q;
                    for (int rep = 0; rep < 3; rep++) {
                        cpam::cpam_query_nodes_touched = 0;
                        h = 0;
                        auto res = CPAMBB::knn(cpambb_master, current_knn[i], K);
                        cnt = res.size();
                        while(!res.empty()) { h += res.top().first.id; res.pop(); }
                    }
                    double q_ms = (t_q.stop() * 1000.0) / 3.0;
                    auto mems = mem_cpambb(cpambb_master);
                    size_t est_nodes = cpam::cpam_query_nodes_touched.load();
                    fout << "CPAMBB | " << qk << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << est_nodes << "\n";
                }
            }
            if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                for(int j=0; j<6; j++) {
                    if (run_algo == names[j] || run_algo == "all") {
                        for(size_t i=0; i<100; i++) {
                            size_t cnt = 0;
                            size_t h = 0;
                            cpam::timer t_q;
                            for (int rep = 0; rep < 3; rep++) {
                                h = 0;
                                auto res = rlog_master[j]->knn_report(current_knn[i], K);
                                cnt = res.size();
                                for(const auto& v : res) h += v.second;
                            }
                            double q_ms = (t_q.stop() * 1000.0) / 3.0;
                            size_t est_nodes = rlog_master[j]->query_lookup_count;
                            auto mems = mem_rlog(*rlog_master[j], rlog_global_history[j]);
                            fout << names[j] << " | " << qk << " | " << i << " | " << cnt << " | " << h << " | " << q_ms << " | " << mems.first << " | " << mems.second << " | " << est_nodes << "\n";
                        }
                    }
                }
            }
        }
        fout << std::flush;
        fout.close();
    };

    size_t master_seq = 0;
    size_t overlap_count = 0;
    size_t commit_count = 0;
    double accumulated_commit_ms = 0.0;
    std::cout << "[Event-Driven] Executing Workload with step size: " << q_step << std::endl;
    
    // Warm up the parallel thread pool to prevent the massive phantom "cold-start" spike on the first commit
    parlay::parallel_for(0, 10000, [](size_t i) { volatile size_t x = i; (void)x; });
    
    ofstream commit_fout("verification_results/CommitLog_" + run_algo + ".txt");
    commit_fout << "Algo | CS_ID | Year | Adds | Rems | Time_ms | Total_Index_MB | Delta_Index_MB | Delta_Data_MB\n";
    
    // 核心诉求：把 Base Build 的巨大耗时，作为 X=0 (2018年初) 的第一个点强行打入日志！
    commit_fout << run_algo << " | 0 | " << start_year << " | " << P_base_set.size() << " | 0 | " << base_build_ms << " | 0 | 0 | 0\n";
    
    CommitStream stream(start_year, dataset_dir, commits_dir, end_year);
    Event ev;
    int ev_year;
    
    while(stream.next_changeset(ev, ev_year)) {
        auto& ops = ev.ops;
        int ev_year_log = ev_year;
        size_t cs_id_log = ev.changeset_id;
            parlay::sequence<geobase::Point> adds, rems;
            vector<Value> boost_adds;
            
            // Collapse operations per node within the changeset to avoid duplicate keys in a batch
            std::unordered_map<size_t, geobase::Point> initial_state;
            std::unordered_map<size_t, Operation> final_ops;
            
            for (auto& op : ops) {
                if (initial_state.find(op.node_id) == initial_state.end() && active_nodes.find(op.node_id) != active_nodes.end()) {
                    initial_state[op.node_id] = active_nodes[op.node_id];
                }
                final_ops[op.node_id] = op;
            }
            
            for (auto& kv : final_ops) {
                auto& op = kv.second;
                bool existed = initial_state.find(op.node_id) != initial_state.end();
                
                if (op.op_type == 'D') {
                    if (existed) {
                        rems.push_back(initial_state[op.node_id]);
                        active_nodes.erase(op.node_id);
                    }
                } else { // I or U
                    if (existed) {
                        auto old_pt = initial_state[op.node_id];
                        if (std::abs(old_pt.x - op.lon) < 0.001 && std::abs(old_pt.y - op.lat) < 0.001) {
                            continue; // Skip NO-OP spatial updates!
                        }
                        rems.push_back(old_pt);
                    }
                    adds.push_back(geobase::Point(op.node_id, op.lon, op.lat));
                    active_nodes[op.node_id] = geobase::Point(op.node_id, op.lon, op.lat);
                    all_ever_existed_nodes.push_back(geobase::Point(op.node_id, op.lon, op.lat));
                    if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                        boost_adds.push_back(make_pair(BoostPoint(op.lon, op.lat), op.node_id));
                    }
                }
            }
            
            double delta_data_mb = (adds.size() * 24.0 - rems.size() * 24.0) / (1024.0 * 1024.0);

            if (run_algo == "MVZD" || run_algo == "all") {
                auto prev_mems = mem_mvzd(); double prev_mem = prev_mems.first + prev_mems.second;
                cpam::timer t_c;
                auto sorted_adds = geobase::get_sorted_points(adds);
                auto sorted_rems = geobase::get_sorted_points(rems);
                mvzd_master = mvzd_tree->commit(mvzd_master, sorted_adds, sorted_rems);
                double cur_ms = t_c.stop() * 1000.0;
                accumulated_commit_ms += cur_ms;
                mvzd_global_history.push_back(mvzd_master);
                double new_mem = mem_mvzd().first + mem_mvzd().second;
                commit_fout << "MVZD | " << cs_id_log << " | " << ev_year_log << " | " << adds.size() << " | " << rems.size() << " | " << cur_ms << " | " << new_mem << " | " << (new_mem - prev_mem) << " | " << delta_data_mb << "\n";
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                auto prev_mems = mem_cpambb(cpambb_master); double prev_mem = prev_mems.first + prev_mems.second;
                cpam::timer t_c;
                auto v_cpambb = CPAMBB::map_insert(adds, cpambb_master);
                if (rems.size() > 0) v_cpambb = CPAMBB::map_delete(rems, v_cpambb);
                double cur_ms = t_c.stop() * 1000.0;
                accumulated_commit_ms += cur_ms;
                cpambb_master = v_cpambb;
                double new_mem = mem_cpambb(cpambb_master).first + mem_cpambb(cpambb_master).second;
                commit_fout << "CPAMBB | " << cs_id_log << " | " << ev_year_log << " | " << adds.size() << " | " << rems.size() << " | " << cur_ms << " | " << new_mem << " | " << (new_mem - prev_mem) << " | " << delta_data_mb << "\n";
            }
            if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                for(int j=0; j<6; j++) {
                    if (run_algo == names[j] || run_algo == "all") {
                        auto prev_mems = mem_rlog(*rlog_master[j], rlog_global_history[j]); double prev_mem = prev_mems.first + prev_mems.second;
                        cpam::timer t_c;
                        RlogBranch v_rlog; 
                        v_rlog.base_snapshot = rlog_master[j]->snapshot;
                        v_rlog.insert_log = boost_adds;
                        
                        // Map deletions to BoostPoint
                        for (auto& r_pt : rems) {
                            v_rlog.remove_log.push_back(make_pair(BoostPoint(r_pt.x, r_pt.y), r_pt.id));
                        }
                        
                        rlog_master[j]->merge(v_rlog);
                        rlog_master[j]->check_and_compact(ev_year_log);
                        
                        double cur_ms = t_c.stop() * 1000.0;
                        accumulated_commit_ms += cur_ms;
                        rlog_global_history[j].push_back(v_rlog);
                        double new_mem = mem_rlog(*rlog_master[j], rlog_global_history[j]).first + mem_rlog(*rlog_master[j], rlog_global_history[j]).second;
                        commit_fout << names[j] << " | " << cs_id_log << " | " << ev_year_log << " | " << adds.size() << " | " << rems.size() << " | " << cur_ms << " | " << new_mem << " | " << (new_mem - prev_mem) << " | " << delta_data_mb << "\n";
                    }
                }
            }

            master_seq++;
            commit_count++;
            
            if (commit_count > 0 && commit_count % q_step == 0) {
                std::cout << "[Event-Driven] Executing Checkpoint " << commit_count << "..." << std::endl;
                string fname = "verification_results/Checkpoint_" + to_string(commit_count) + "_" + run_algo + ".txt";
                
                ofstream fout(fname);
                fout << "Algo | QType | QID | Count | Hash | Time_ms | Mem_Tree_MB | Mem_Delta_MB | Nodes\n";
                fout.close();
                
                run_hybrid_queries(fname);
                
                accumulated_commit_ms = 0.0;
            }
    }
    
    ofstream end_flag("verification_results/2026_EndQuery_" + run_algo + ".txt");
    end_flag << "Overlaps: " << overlap_count << std::endl;
    return 0;
}
