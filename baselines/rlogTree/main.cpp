#include<bits/stdc++.h>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/index/detail/rtree/utilities/view.hpp>
#include "geo/point.hpp"
#include "geo/operations.hpp"
#include "geo/io.hpp"
#include "cpam/parse_command_line.h"
#include "helper/time_loop.h"

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> Point;
typedef pair<Point, size_t> Value; // (point, ID)
typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;
size_t maxSize = 100;

template<typename PT>
auto convert_P(PT &P){
    vector<Value> ret(P.size());
    parlay::parallel_for(0, P.size(), [&](size_t i){
        ret[i] = make_pair(Point(P[i].x, P[i].y), P[i].id);
    });
    return ret;
}

// -------------------------------------------------------------------
// Rlog Tree Wrapper
// Implements Snapshot + Append-only Delta Log
// -------------------------------------------------------------------
class RlogTree {
public:
    shared_ptr<RTree> snapshot;
    vector<Value> insert_log;
    vector<Value> remove_log;
    int version_depth;
    int threshold;

    RlogTree(int t = 4) : version_depth(0), threshold(t) {
        snapshot = make_shared<RTree>();
    }

    // Build initial base tree
    void build_base(const vector<Value>& pts) {
        snapshot = make_shared<RTree>(pts.begin(), pts.end());
        insert_log.clear();
        remove_log.clear();
        version_depth = 0;
    }

    // Copy Constructor: O(1) Fork! (just shared_ptr ref count + vector copies)
    RlogTree(const RlogTree& other) {
        snapshot = other.snapshot;
        insert_log = other.insert_log;
        remove_log = other.remove_log;
        version_depth = other.version_depth;
        threshold = other.threshold;
    }

    RlogTree& operator=(const RlogTree& other) {
        if (this != &other) {
            snapshot = other.snapshot;
            insert_log = other.insert_log;
            remove_log = other.remove_log;
            version_depth = other.version_depth;
            threshold = other.threshold;
        }
        return *this;
    }

    void commit_inserts(const vector<Value>& pts) {
        insert_log.insert(insert_log.end(), pts.begin(), pts.end());
    }

    void commit_removes(const vector<Value>& pts) {
        remove_log.insert(remove_log.end(), pts.begin(), pts.end());
    }

    void compact() {
        // 1. Deep copy the snapshot (EXPENSIVE)
        auto new_snap = make_shared<RTree>(*snapshot);
        // 2. Apply logs
        for(auto& p : remove_log) new_snap->remove(p);
        new_snap->insert(insert_log.begin(), insert_log.end());
        // 3. Reset
        snapshot = new_snap;
        insert_log.clear();
        remove_log.clear();
        version_depth = 0;
    }

    void merge(const RlogTree& branch, size_t original_log_size) {
        size_t new_inserts_size = branch.insert_log.size() - original_log_size;
        if (new_inserts_size > 0) {
            insert_log.insert(
                insert_log.end(), 
                branch.insert_log.begin() + original_log_size, 
                branch.insert_log.end()
            );
        }
        version_depth++;

        if (version_depth >= threshold) {
            compact();
        }
    }

    vector<Value> range_report(const geobase::Bounding_Box& q) const {
        vector<Value> result;
        bg::model::box<Point> box(Point(q.first.x, q.first.y), Point(q.second.x, q.second.y));
        snapshot->query(bgi::intersects(box), std::back_inserter(result));
        
        std::set<size_t> removed_ids;
        for (const auto& val : remove_log) {
            if (bg::within(val.first, box)) {
                 removed_ids.insert(val.second);
            }
        }
        
        vector<Value> final_result;
        for (const auto& val : result) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                final_result.push_back(val);
            }
        }
        
        for (const auto& val : insert_log) {
            if (bg::within(val.first, box)) {
                final_result.push_back(val);
            }
        }
        return final_result;
    }
};



#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>

void run_evolution_test(string dir, string q0, string q1, string q2, int threshold) {
    auto load_q = [](string f) {
        parlay::sequence<geobase::Bounding_Box> queries;
        if (f != "") {
            auto [cnt, q] = geobase::read_range_query(f, 8, maxSize);
            queries = q.substr(0, std::min((size_t)100, q.size()));
        }
        return queries;
    };
    auto queries_S = load_q(q0);
    auto queries_M = load_q(q1);
    auto queries_L = load_q(q2);

    string base_file = dir + "/base_upto_2020.txt";
    ifstream fin_base(base_file);
    if (!fin_base.is_open()) return;
    parlay::sequence<geobase::Point> P_base;
    size_t tid; double tx, ty;
    while (fin_base >> tid >> tx >> ty) {
        if (tx < 0 || ty < 0) continue;
        P_base.emplace_back(geobase::Point(tid, tx*1000000, ty*1000000));
    }
    auto P_base_conv = convert_P(P_base);
    
    RlogTree master_ver(threshold);
    master_ver.build_base(P_base_conv);
    cout << "[INFO] Base tree loaded. Nodes: " << P_base.size() << endl;

    auto run_queries = [&](const RlogTree& tree, const parlay::sequence<geobase::Bounding_Box>& queries, string name) {
        if (queries.size() == 0) return;
        ofstream fout("verification_results/raw_" + name + ".txt");
        fout << "Tree Size: " << tree.snapshot->size() + tree.insert_log.size() - tree.remove_log.size() << "
";
        for(size_t j=0; j<queries.size(); j++) {
            parlay::internal::timer t; t.start();
            auto out = tree.range_report(queries[j]);
            double elapsed_ms = t.stop() * 1000.0;
            fout << "Q" << j << ": " << out.size() << " " << fixed << setprecision(6) << elapsed_ms << "
";
        }
        fout.close();
    };

    for (int year = 2021; year <= 2025; year++) {
        vector<string> uids = {"66497", "3442203", "352700", "others"};
        if (year == 2022) uids = {"66497", "822681", "670875", "others"};
        if (year == 2023) uids = {"66497", "670875", "14496019", "others"};
        if (year == 2024) uids = {"670875", "6518572", "14570006", "others"};
        if (year == 2025) uids = {"1983103", "437598", "9560399", "others"};
        
        RlogTree merged_ver = master_ver;
        for (size_t i = 0; i < uids.size(); i++) {
            string branch_file = dir + "/" + to_string(year) + "_branch_" + to_string(i+1) + (uids[i] == "others" ? "_others.txt" : "_uid_" + uids[i] + ".txt");
            ifstream fin_branch(branch_file);
            parlay::sequence<geobase::Point> P_branch;
            while (fin_branch >> tid >> tx >> ty) {
                if (tx < 0 || ty < 0) continue;
                P_branch.emplace_back(geobase::Point(tid, tx*1000000, ty*1000000));
            }
            auto P_branch_conv = convert_P(P_branch);

            RlogTree branch_ver = master_ver;
            branch_ver.commit_inserts(P_branch_conv);
            merged_ver.merge(branch_ver, master_ver.insert_log.size());
            
            string pfx = "rlog" + to_string(threshold) + "_" + to_string(year) + "_B" + to_string(i+1);
            run_queries(merged_ver, queries_S, pfx + "_S");
            run_queries(merged_ver, queries_M, pfx + "_M");
            run_queries(merged_ver, queries_L, pfx + "_L");
        }
        master_ver = merged_ver;
        cout << "Year " << year << " completed for Rlog " << threshold << endl;
    }
}

void run(int argc, char** argv){
    cpam::commandLine cmd(argc, argv, "[-i <Path-to-Input>] [-o <Path-to-Output>] [-t <Task-Name>] [-a <Algorithm-Name>] "
									  "[-b <Path-to-Batch-file>] [-bf <batch-fraction>] "
									  "[-r <Path-to-Range-Query>] [-real <Is-Real-Dataset?>] "
									  "[-mv <Dir-to-Multi-Version>] [-s <Single-Point-Query>]"
									  );
    string task = cmd.getOptionValue("-t");
    
    int threshold = cmd.getOptionIntValue("-th", 4);
    if (task == "bhutan-merge"){

        string dir = cmd.getOptionValue("-mv");
        run_evolution_test(dir, cmd.getOptionValue("-r0", ""), cmd.getOptionValue("-r1", ""), cmd.getOptionValue("-r2", ""), threshold);
        return;
    }
}

int main(int argc, char** argv){
    run(argc, argv);
    return 0;
}
