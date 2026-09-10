#include<bits/stdc++.h>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/index/detail/rtree/utilities/view.hpp>
#include "geo/point.hpp"
#include "geo/operations.hpp"
#include "geo/io.hpp"
#include "cpam/parse_command_line.h"
#include "helper/time_loop.h"

#define TEST

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define point and R-tree value type
typedef bg::model::point<double, 2, bg::cs::cartesian> Point;
typedef pair<Point, size_t> Value; // (point, ID)
typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;

geobase::Bounding_Box largest_mbr;  //  largest mbr 
size_t maxSize = 100;   //  max points in a region

/* PT is a sequence of points */
template<typename PT>
auto convert_P(PT &P){
    vector<Value> ret(P.size());
    parlay::parallel_for(0, P.size(), [&](size_t i){
        ret[i] = make_pair(Point(P[i].x, P[i].y), P[i].id);
    });
    return ret;
}

/* T is a sequence of points */
template<typename T>
auto convert_Q(T &Q){
    vector<pair<Point, Point> > ret(Q.size());
    parlay::parallel_for(0, Q.size(), [&](size_t i){
        auto min_pt = Q[i].first;
        auto max_pt = Q[i].second;
        ret[i] = make_pair(Point(min_pt.x, min_pt.y), Point(max_pt.x, max_pt.y));
    });
    return ret;
}

// Get resident set size (RSS) in MB on Linux
double getProcessRSS_MB() {
    ifstream stat_stream("/proc/self/statm");
    size_t size = 0, resident = 0;
    stat_stream >> size >> resident;
    stat_stream.close();

    long page_size = sysconf(_SC_PAGESIZE); // bytes per page
    size_t rss_bytes = resident * page_size;

    return static_cast<double>(rss_bytes) / (1024.0 * 1024.0); // convert to MB
}

template <typename T, typename DIFF>
auto diff_by_id(T &a, T &b, DIFF &ret_diff){
    // unordered_map<size_t, Point> map_a, map_b;

    // for (const auto &[pt, id]: a){
    //     map_a[id] = pt;
    // }
    // for (const auto &[pt, id]: b){
    //     map_b[id] = pt;
    // }
    
    // for (const auto &[id, pa]: map_a){
    //     auto it = map_b.find(id);
    //     if (it == map_b.end() || !bg::equals(it->second, pa)){
    //         auto cur_p = geobase::Point(id, bg::get<0>(pa), bg::get<1>(pa));
    //         ret_diff.remove_point(cur_p);
    //     }
    // }

    // for (const auto &[id, pb]: map_b){
    //     auto it = map_a.find(id);
    //     if (it == map_a.end() || !bg::equals(it->second, pb)){
    //         auto cur_p = geobase::Point(id, bg::get<0>(pb), bg::get<1>(pb));
    //         ret_diff.add_point(cur_p);
    //     }
    // }

    size_t i = 0, j = 0;

    while (i < a.size() && j < b.size()){
        if (a[i].second == b[j].second){ // same point, check coordinate
            if (bg::equals(a[i].first, b[j].first)){ //  same coordinate, do nothing
                i++, j++;
            }
            else{
                auto cur_pa = geobase::Point(a[i].second, bg::get<0>(a[i].first), bg::get<1>(a[i].first));
                auto cur_pb = geobase::Point(b[j].second, bg::get<0>(b[j].first), bg::get<1>(b[j].first));
                i++, j++;
                ret_diff.remove_point(cur_pa);
                ret_diff.add_point(cur_pb);
            }
        }
        else if (a[i].second < b[j].second){ // first smaller, in A not in B
            auto cur_pa = geobase::Point(a[i].second, bg::get<0>(a[i].first), bg::get<1>(a[i].first));
            ret_diff.remove_point(cur_pa);
            i++;
        }
        else{ //  second smaller, in B not in A
            auto cur_pb = geobase::Point(b[j].second, bg::get<0>(b[j].first), bg::get<1>(b[j].first));
            ret_diff.add_point(cur_pb);
            j++;
        }
    }
    while (i < a.size()){
        auto cur_pa = geobase::Point(a[i].second, bg::get<0>(a[i].first), bg::get<1>(a[i].first));
        ret_diff.add_point(cur_pa);
        i++;
    }
    while (j < b.size()){
        auto cur_pb = geobase::Point(b[j].second, bg::get<0>(b[j].first), bg::get<1>(b[j].first));
        ret_diff.add_point(cur_pb);
        j++;
    }
}

void print_values(const std::vector<Value>& values) {
    for (const auto& val : values) {
        double x = bg::get<0>(val.first);
        double y = bg::get<1>(val.first);
        size_t id = val.second;
        cout << fixed << setprecision(6) << "ID: " << id << " at (" << x << ", " << y << ")\n";
    }
}

template<typename Tree, typename PSet>
void knn_test(Tree &tree, PSet &P, size_t k = 10, size_t q_num = 50000, size_t d = 2){
    /* build tree */
    // q_num = min(q_num, P.size());
    std::vector<Value> res;
    res.reserve(k);

    auto avg_time = time_loop(
        3, 1.0, 
        [&]() {
            res.resize(0);
        },
        [&]() {    
            for (size_t i = 0; i < q_num; i++){
                tree.query(bgi::nearest(P[i].first, k), std::back_inserter(res));
            }
        },
        [&](){
        }
    );

    cout << fixed << setprecision(6) << "[boostRtree-KNN]: " << avg_time << endl;
}

template<typename PT, typename RQ>
void range_report_test(PT P, RQ Q, parlay::sequence<size_t> &cnt){
    auto P_conv = convert_P(P);
    auto Q_conv = convert_Q(Q);

    RTree rtree(P_conv.begin(), P_conv.end());  //  original version 

    parlay::sequence<size_t> rangeReportCnt(Q.size());
    // parlay::sequence<parlay::sequence<Point> > rangeReport(querys.size());
    vector<vector<Value> > rangeReport(Q.size());
    for (size_t i = 0; i < Q.size(); i++){
        // rangeReport[i].resize(cnt[i]);
        rangeReport[i].reserve(cnt[i]);
    }

    for (size_t i = 0; i < Q_conv.size(); i++){
        // print_mbr(querys[i]);
        auto avg_time = time_loop(
            3, 1.0, 
            [&]() {
                rangeReport[i].resize(0);
            },
            [&]() {					
                rtree.query(bgi::intersects(bg::model::box<Point>(Q_conv[i].first, Q_conv[i].second)), std::back_inserter(rangeReport[i]));
                rangeReportCnt[i] = rangeReport[i].size();
            },
            [&](){} );
        if (rangeReportCnt[i] != cnt[i]){
            cout << "[ERROR] Incorrect" << endl;
            cout << rangeReportCnt[i] << " " << cnt[i] << endl;
        }
        else{
            cout << fixed << setprecision(6) << rangeReportCnt[i] << " " << avg_time << endl;
        }
    }
}

template<typename PT, typename RQ>
void spatial_diff_test(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
    auto P_conv = convert_P(P);
    RTree rtree(P_conv.begin(), P_conv.end());  //  original version
    cout << "original version bulkload finished." << endl;
    // cout << rtree.depth() << endl;
    auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
    if (P.size() < max_batch_size) max_batch_size = P.size();

    /* get insert, delete points */
    auto P_test = geobase::shuffle_point(P, max_batch_size);
    cout << "shuffle finished." << endl;
    auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());
    cout << "split finished." << endl;
    auto Q_conv = convert_Q(range_queries);
    cout << "query convert finished." << endl;
    
    for (auto &batch_size: batch_sizes){
        if (batch_size > P.size()) batch_size = P.size();
        cout << "[INFO] Batch Size: " << batch_size << endl;
        auto insert_num = batch_size / 10 * insert_ratio;
        auto delete_num = batch_size / 10 * (10 - insert_ratio);
        
        auto P_insert = P_insert_set.substr(0, insert_num);
        auto P_delete = P_delete_set.substr(0, delete_num);

        /* apply delete, insert */
        auto P_newver = geobase::collect_newver_point(P, P_insert, P_delete);
        auto P_newver_conv = convert_P(P_newver);
        RTree rtree2(P_newver_conv.begin(), P_newver_conv.end());
        cout << "new version bulkload finished." << endl;

        parlay::sequence<size_t> addCnt(range_queries.size());
        parlay::sequence<size_t> removeCnt(range_queries.size());
        // vector<Value> res1, res2;
        parlay::sequence<Value> res1, res2;
        geobase::diff_type ret_diff(maxSize, maxSize);
        res1.reserve(2 * maxSize);
        res2.reserve(2 * maxSize);

        for (size_t i = 0; i < Q_conv.size(); i++){
            auto avg_time = time_loop(
                3, 1.0,
                [&](){},
                [&](){
                        ret_diff.reset(maxSize, maxSize);
                        res1.resize(0);
                        res2.resize(0);
                        rtree.query(bgi::intersects(bg::model::box<Point>(Q_conv[i].first, Q_conv[i].second)), std::back_inserter(res1));
                        rtree2.query(bgi::intersects(bg::model::box<Point>(Q_conv[i].first, Q_conv[i].second)), std::back_inserter(res2));

                        auto id_cmp = [&](auto &lhs, auto &rhs) { return lhs.second < rhs.second; };

                        res1 = parlay::sort(res1, id_cmp);
                        res2 = parlay::sort(res2, id_cmp);

                        diff_by_id(res1, res2, ret_diff);
                        ret_diff.compact();
                        // auto [add, remove] = diff_by_id(res1, res2);
                        addCnt[i] = ret_diff.add.size();
                        removeCnt[i] = ret_diff.remove.size();
                },
                [&]{}
            );
            cout << fixed << setprecision(6) << i << " " << avg_time << endl;
        }

        // auto avg_time = time_loop(
        //     3, 1.0,
        //     [&](){},
        //     [&](){
        //         for (size_t i = 0; i < Q_conv.size(); i++){
        //             ret_diff.reset(maxSize, maxSize);
        //             res1.resize(0);
        //             res2.resize(0);
        //             rtree.query(bgi::intersects(bg::model::box<Point>(Q_conv[i].first, Q_conv[i].second)), std::back_inserter(res1));
        //             rtree2.query(bgi::intersects(bg::model::box<Point>(Q_conv[i].first, Q_conv[i].second)), std::back_inserter(res2));

        //             auto id_cmp = [&](auto &lhs, auto &rhs) { return lhs.second < rhs.second; };

        //             res1 = parlay::sort(res1, id_cmp);
        //             res2 = parlay::sort(res2, id_cmp);

        //             diff_by_id(res1, res2, ret_diff);
        //             ret_diff.compact();
        //             // auto [add, remove] = diff_by_id(res1, res2);
        //             addCnt[i] = ret_diff.add.size();
        //             removeCnt[i] = ret_diff.remove.size();
        //         }
        //     },
        //     [&]{}
        // );
        // cout << fixed << setprecision(6) << "[Boost-Rtree] spatial-diff time (avg): " << avg_time << endl;

        #ifdef TEST
        string file_name = "output/rtree_spatial_diff-" + to_string(batch_size); 
        ofstream spatialDiffOut(file_name);
        for (size_t i = 0; i < range_queries.size(); i++){
            spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
        }
        #endif
    }
}

template<typename PT>
void multi_version_test(PT P, string dir, int start_year = 14, int version_num = 5){
    auto cur_year = start_year;
    
        vector<parlay::sequence<geobase::Point>> P_delete(version_num), P_insert(version_num), P_update(version_num), P_updove(version_num);
    vector<vector<Value>> P_delete_conv(version_num), P_insert_conv(version_num);
    

    for (auto i = 0; i != version_num; i++){
        auto delete_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-delete.txt";
        auto insert_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-insert.txt";
        auto update_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update.txt";
        auto updove_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update_remove.txt";

        ifstream fin_delete(delete_file_name);
        ifstream fin_insert(insert_file_name);
        ifstream fin_update(update_file_name);
        ifstream fin_updove(updove_file_name);

        auto delete_mbr = read_pts(P_delete[i], fin_delete, 1);
        auto insert_mbr = read_pts(P_insert[i], fin_insert, 1);
        auto update_mbr = read_pts(P_update[i], fin_update, 1);
        auto updove_mbr = read_pts(P_updove[i], fin_updove, 1);

        P_delete[i].append(P_updove[i]);
        P_insert[i].append(P_update[i]);

        P_delete_conv[i] = convert_P(P_delete[i]);
        P_insert_conv[i] = convert_P(P_insert[i]);

        cur_year += 1;
        delete_mbr = insert_mbr; // useless, just remove warning
        update_mbr = updove_mbr;
    }

    auto P_conv = convert_P(P); // transfer to boost format

    /* Build initial version */
    auto build_avg = time_loop(
        3, 1.0, [&]() {},
        [&]() {
            RTree tst_rtree(P_conv.begin(), P_conv.end());
        },
        [&](){});

    vector<RTree> rtrees(version_num + 1);

    auto before_mem = getProcessRSS_MB();
    rtrees[0].insert(P_conv.begin(), P_conv.end());
    auto after_mem = getProcessRSS_MB();

    cout << "[Rtree init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;
    cout << "[Rtree init build memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;

    for (auto i = 0; i < version_num; i++){
        cout << "dealing with version " << i + 1 << ":" << endl;

        auto commit_avg = time_loop(
            3, 1.0, 
            [&]() {
                rtrees[i + 1].clear();
            },
            [&]() {
                before_mem = getProcessRSS_MB();
                rtrees[i + 1] = rtrees[i];
                for (auto &pt: P_delete_conv[i]){
                    rtrees[i + 1].remove(pt);
                }
                rtrees[i + 1].insert(P_insert_conv[i].begin(), P_insert_conv[i].end());
                after_mem = getProcessRSS_MB();
            },
            [&](){
            });
        cout << rtrees[i + 1].size() << endl;
        cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;
        cout << "[Rtree new version memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;
    }
}




#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>


#include <sys/stat.h>

void run_evolution_test(string dir, string q0, string q1, string q2) {
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
    RTree master_ver(P_base_conv.begin(), P_base_conv.end());
    cout << "[INFO] Base tree loaded. Nodes: " << P_base.size() << endl;

    auto run_queries = [&](const RTree& tree, const parlay::sequence<geobase::Bounding_Box>& queries, string name) {
        if (queries.size() == 0) return;
        ofstream fout("verification_results/raw_" + name + ".txt");
        fout << "Tree Size: " << tree.size() << "
";
        for(size_t j=0; j<queries.size(); j++) {
            bg::model::box<Point> box(Point(queries[j].first.x, queries[j].first.y), Point(queries[j].second.x, queries[j].second.y));
            std::vector<Value> result;
            parlay::internal::timer t; t.start();
            tree.query(bgi::intersects(box), std::back_inserter(result));
            double elapsed_ms = t.stop() * 1000.0;
            fout << "Q" << j << ": " << result.size() << " " << fixed << setprecision(6) << elapsed_ms << "
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
        
        RTree merged_ver = master_ver;
        for (size_t i = 0; i < uids.size(); i++) {
            string branch_file = dir + "/" + to_string(year) + "_branch_" + to_string(i+1) + (uids[i] == "others" ? "_others.txt" : "_uid_" + uids[i] + ".txt");
            ifstream fin_branch(branch_file);
            parlay::sequence<geobase::Point> P_branch;
            while (fin_branch >> tid >> tx >> ty) {
                if (tx < 0 || ty < 0) continue;
                P_branch.emplace_back(geobase::Point(tid, tx*1000000, ty*1000000));
            }
            auto P_branch_conv = convert_P(P_branch);

            RTree branch_ver = master_ver;
            branch_ver.insert(P_branch_conv.begin(), P_branch_conv.end());
            
            std::vector<Value> A(master_ver.begin(), master_ver.end());
            std::vector<Value> B(branch_ver.begin(), branch_ver.end());
            auto cmp = [](const Value& a, const Value& b) { return a.second < b.second; };
            std::sort(A.begin(), A.end(), cmp);
            std::sort(B.begin(), B.end(), cmp);
            std::vector<Value> added, removed;
            std::set_difference(B.begin(), B.end(), A.begin(), A.end(), std::back_inserter(added), cmp);
            std::set_difference(A.begin(), A.end(), B.begin(), B.end(), std::back_inserter(removed), cmp);
            for(auto& v : removed) merged_ver.remove(v);
            merged_ver.insert(added.begin(), added.end());
            
            string pfx = "boost_" + to_string(year) + "_B" + to_string(i+1);
            run_queries(merged_ver, queries_S, pfx + "_S");
            run_queries(merged_ver, queries_M, pfx + "_M");
            run_queries(merged_ver, queries_L, pfx + "_L");
        }
        master_ver = merged_ver;
        cout << "Year " << year << " completed for Boost" << endl;
    }
}

void run(int argc, char** argv){
    cpam::commandLine cmd(argc, argv, "[-i <Path-to-Input>] [-o <Path-to-Output>] [-t <Task-Name>] [-a <Algorithm-Name>] "
									  "[-b <Path-to-Batch-file>] [-bf <batch-fraction>] "
									  "[-r <Path-to-Range-Query>] [-real <Is-Real-Dataset?>] "
									  "[-mv <Dir-to-Multi-Version>] [-s <Single-Point-Query>]"
									  );
	if (!cmd.getOption("-t")){ 
        cout << "[ERROR]: <Task-Name> is not specified." << endl;
		return;
	}
	if (!cmd.getOption("-i")){
		cout << "[ERROR]: <Path-to-Input> is not specified." << endl;
		return;
	}

	string task = cmd.getOptionValue("-t");
	string input_file = cmd.getOptionValue("-i");
	int is_real = cmd.getOptionIntValue("-real", 0);

	/* read input file */
	ifstream fin(input_file);
	parlay::sequence<geobase::Point> P;
	geobase::read_pts(P, fin, is_real);	//	change to true if id is contained.

    // auto P_conv = convert_P(P); // get point format for boost
    // RTree rtree(P_conv.begin(), P_conv.end());
    // cout << "bulk-load finished." << endl;

    auto P_conv = convert_P(P); // transfer to boost format
    RTree rtree(P_conv.begin(), P_conv.end());  //  original version

    if (task == "knn"){
        parlay::sequence<size_t> k_vals = {
            1,
            2,
            5,
            10,
            20,
            50,
            100
        };
        for (auto &k: k_vals){
            cout << "[INFO] k = " << k << endl;
            knn_test(rtree, P_conv, k);
        }    
    }
    
    if (task == "range-report"){
		string query_file = cmd.getOptionValue("-r");
		auto [cnt, queries] = geobase::read_range_query(query_file, 8, maxSize);
        cout << "boost-Rtree: " << endl;
        range_report_test(P, queries, cnt);
    }

    if (task == "spatial-diff"){
		string query_file = cmd.getOptionValue("-r");
		auto [cnt, queries] = geobase::read_range_query(query_file, 8, maxSize);
		// size_t batch_size = 10000;
        // queries = queries.substr(0, 1);


		parlay::sequence<size_t> batch_sizes = {
            // 10000,
            // 20000,
            // 50000,
            // 100000,
            // 200000,
            // 500000,
            // 1000000,
            // 2000000,
            // 5000000,
            10000000,
            // 20000000, 
            // 50000000,
            // 100000000
        };
        // parlay::sequence<size_t> ratios = {1, 2, 3, 4, 5, 6, 7, 8, 9};
		// queries = queries.substr(6148, 1);
        size_t ratio = 5;
        cout << "[BoostRtree]" << endl;
        spatial_diff_test(P, queries, batch_sizes, ratio);
	}

    	/* multi-version test */
	
    if (task == "bhutan-merge"){
        string dir = cmd.getOptionValue("-mv");
        run_evolution_test(dir, cmd.getOptionValue("-r0", ""), cmd.getOptionValue("-r1", ""), cmd.getOptionValue("-r2", ""));
        return;
    }
    if (task == "multi-version-test"){

		if (!cmd.getOption("-mv")){
			cout << "[ERROR]: <Dir-to-Multi-Version> is not specified." << endl;
		}
		else{
            string mv_dir = cmd.getOptionValue("-mv");
            // multi_version_test(P, mv_dir, 18, 7);
            multi_version_test(P, mv_dir, 14, 11);   // Japan
		}
		return;
	} 

}

int main(int argc, char** argv){
    run(argc, argv);
    return 0;
}