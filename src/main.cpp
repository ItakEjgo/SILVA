#include <iostream>
#include <fstream>
#include <string>

#include <cpam/parse_command_line.h>
#include <silva/geo/io.hpp>
#include <silva/core/global_config.hpp>

// Include the test suites
#include "../tests/test_rlog.hpp"
#include "../tests/test_pkd.hpp"
#include "../tests/test_pkdlog.hpp"
#include "../tests/test_boost.hpp"
#include "../tests/test_knn_correctness.hpp"
#include "../tests/test_mvq.hpp"
#include "../tests/test_pacz.hpp"

using namespace std;

void line_splitter() {
    cout << "-------------------------------------------------------" << endl;
}

void run(int argc, char** argv) {
    cpam::commandLine cmd(argc, argv, "[-i <Path-to-Input>] [-u <Path-to-Update>] [-t <Task-Name>] [-a <Algorithm-Name>] "
                                      "[-r <Path-to-Range-Query>] [-real <Is-Real-Dataset?>] "
                                      "[-br <Batch-Ratios>] [-k <KNN-K>] [-qn <Query-Num>]");
    
    if (!cmd.getOption("-t") || !cmd.getOption("-i") || !cmd.getOption("-a")) {
        cout << "[ERROR]: Missing required arguments: -t <Task-Name>, -i <Path-to-Input>, -a <Algorithm-Name>" << endl;
        return;
    }

    string task = cmd.getOptionValue("-t");
    string algo = cmd.getOptionValue("-a");
    string input_file = cmd.getOptionValue("-i");
    string update_file = cmd.getOptionValue("-u", "");
    int is_real = cmd.getOptionIntValue("-real", 0);

    // Read input file
    ifstream fin(input_file);
    if (!fin.is_open()) {
        cout << "[ERROR]: Cannot open input file: " << input_file << endl;
        return;
    }

    parlay::sequence<geobase::Point> P_base;
    mvq::Config::get().largest_mbr = geobase::read_pts(P_base, fin, is_real);
    
    parlay::sequence<geobase::Point> P_update;
    if (task == "batch-insert" && update_file != "") {
        ifstream fin_u(update_file);
        if (!fin_u.is_open()) {
            cout << "[ERROR]: Cannot open update file: " << update_file << endl;
            return;
        }
        geobase::read_pts(P_update, fin_u, is_real);
    }

    if (task == "debug") {
        cout << "total points base: " << P_base.size() << endl;
        if (update_file != "") cout << "total points update: " << P_update.size() << endl;
        return;
    }

    // Batch ratios definition
    string batch_ratios_str = cmd.getOptionValue("-br", "0.01,0.1,0.25,0.5,1.0");
    double compact_p = cmd.getOptionDoubleValue("-p", 1.0);
    parlay::sequence<double> batch_ratios;
    size_t pos = 0;
    while ((pos = batch_ratios_str.find(",")) != string::npos) {
        batch_ratios.push_back(stod(batch_ratios_str.substr(0, pos)));
        batch_ratios_str.erase(0, pos + 1);
    }
    batch_ratios.push_back(stod(batch_ratios_str));

    if (task == "build") {
        if (algo == "mvq" || algo == "combined") ZDTest::build_test(P_base);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::build_test(P_base);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::build_test(P_base);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::build_test(P_base);
        if (algo == "combined") line_splitter();
        if (algo == "pkdlog" || algo == "combined") PKDLogTest::build_test(P_base);
        if (algo == "combined") line_splitter();
        if (algo == "boost" || algo == "combined") BoostTest::build_test(P_base);
    } 
    else if (task == "batch-insert") {
        if (update_file == "") {
            cout << "[ERROR]: batch-insert requires -u <Path-to-Update>" << endl;
            return;
        }
        if (algo == "mvq" || algo == "combined") ZDTest::batch_insert_test(P_base, P_update, batch_ratios);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::batch_insert_test(P_base, P_update, batch_ratios);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::batch_insert_test(P_base, P_update, batch_ratios, compact_p);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::batch_insert_test(P_base, P_update, batch_ratios);
        if (algo == "combined") line_splitter();
        if (algo == "pkdlog" || algo == "combined") PKDLogTest::batch_insert_test(P_base, P_update, batch_ratios, compact_p);
        if (algo == "combined") line_splitter();
        if (algo == "boost" || algo == "combined") BoostTest::batch_insert_test(P_base, P_update, batch_ratios);
    }
    else if (task == "batch-delete") {
        if (algo == "mvq" || algo == "combined") ZDTest::batch_delete_test(P_base, batch_ratios);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::batch_delete_test(P_base, batch_ratios);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::batch_delete_test(P_base, batch_ratios, compact_p);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::batch_delete_test(P_base, batch_ratios);
        if (algo == "combined") line_splitter();
        if (algo == "pkdlog" || algo == "combined") PKDLogTest::batch_delete_test(P_base, batch_ratios, compact_p);
        if (algo == "combined") line_splitter();
        if (algo == "boost" || algo == "combined") BoostTest::batch_delete_test(P_base, batch_ratios);
    }
    else if (task == "range-count") {
        string count_qry_file = cmd.getOptionValue("-r", "range_count.qry");
        auto q_tuple = geobase::read_range_query(count_qry_file, 4, mvq::Config::get().maxSize);
        auto querys = std::get<1>(q_tuple);
        auto cnt = std::get<0>(q_tuple);
        if (algo == "mvq" || algo == "combined") ZDTest::range_count_test(P_base, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::range_count_test(P_base, querys, cnt);
    }
    else if (task == "range-report") {
        string report_qry_file = cmd.getOptionValue("-r", "range_report.qry");
        auto q_tuple = geobase::read_range_query(report_qry_file, 8, mvq::Config::get().maxSize);
        auto querys = std::get<1>(q_tuple);
        auto cnt = std::get<0>(q_tuple);
        if (algo == "mvq" || algo == "combined") ZDTest::range_report_test(P_base, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::range_report_test(P_base, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::range_report_test(P_base, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::range_report_test(P_base, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "pkdlog" || algo == "combined") PKDLogTest::range_report_test(P_base, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "boost" || algo == "combined") BoostTest::range_report_test(P_base, querys, cnt);
    }
    else if (task == "knn") {
        size_t k = cmd.getOptionIntValue("-k", 10);
        size_t q_num = cmd.getOptionIntValue("-qn", 50000);
        if (algo == "mvq" || algo == "combined") ZDTest::knn_test(P_base, k, q_num);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::knn_test(P_base, k, q_num);
    }
    else if (task == "knn-verify") {
        KNNVerify::run_verification(P_base);
    }
    else {
        cout << "[ERROR]: Unknown task: " << task << endl;
    }
}

int main(int argc, char** argv) {
    run(argc, argv);
    return 0;
}
