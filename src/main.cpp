#include <iostream>
#include <fstream>
#include <string>

#include <cpam/parse_command_line.h>
#include <silva/geo/io.hpp>
#include <silva/core/global_config.hpp>

// Include the test suites
#include "../tests/test_rlog.hpp"
#include "../tests/test_pkd.hpp"
#include "../tests/test_mvq.hpp"
#include "../tests/test_pacz.hpp"

using namespace std;

void line_splitter() {
    cout << "-------------------------------------------------------" << endl;
}

void run(int argc, char** argv) {
    cpam::commandLine cmd(argc, argv, "[-i <Path-to-Input>] [-t <Task-Name>] [-a <Algorithm-Name>] "
                                      "[-r <Path-to-Range-Query>] [-real <Is-Real-Dataset?>] "
                                      "[-k <KNN-K>] [-qn <Query-Num>]");
    
    if (!cmd.getOption("-t") || !cmd.getOption("-i") || !cmd.getOption("-a")) {
        cout << "[ERROR]: Missing required arguments: -t <Task-Name>, -i <Path-to-Input>, -a <Algorithm-Name>" << endl;
        return;
    }

    string task = cmd.getOptionValue("-t");
    string algo = cmd.getOptionValue("-a");
    string input_file = cmd.getOptionValue("-i");
    int is_real = cmd.getOptionIntValue("-real", 0);

    // Read input file
    ifstream fin(input_file);
    if (!fin.is_open()) {
        cout << "[ERROR]: Cannot open input file: " << input_file << endl;
        return;
    }

    parlay::sequence<geobase::Point> P;
    mvq::Config::get().largest_mbr = geobase::read_pts(P, fin, is_real);

    if (task == "debug") {
        cout << "total points: " << P.size() << endl;
        return;
    }

    // Batch sizes definition
    parlay::sequence<size_t> batch_sizes = {
        10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000,
        10000000, 20000000, 50000000, 100000000
    };

    if (task == "build") {
        if (algo == "mvq" || algo == "combined") ZDTest::build_test(P);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::build_test(P);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::build_test(P);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::build_test(P);
    } 
    else if (task == "batch-insert") {
        if (algo == "mvq" || algo == "combined") ZDTest::batch_insert_test(P, batch_sizes);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::batch_insert_test(P, batch_sizes);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::batch_insert_test(P, batch_sizes);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::batch_insert_test(P, batch_sizes);
    }
    else if (task == "batch-delete") {
        if (algo == "mvq" || algo == "combined") ZDTest::batch_delete_test(P, batch_sizes);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::batch_delete_test(P, batch_sizes);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::batch_delete_test(P, batch_sizes);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::batch_delete_test(P, batch_sizes);
    }
    else if (task == "range-count") {
        string count_qry_file = cmd.getOptionValue("-r", "range_count.qry");
        auto q_tuple = geobase::read_range_query(count_qry_file, 4, mvq::Config::get().maxSize);
        auto querys = std::get<1>(q_tuple);
        auto cnt = std::get<0>(q_tuple);
        if (algo == "mvq" || algo == "combined") ZDTest::range_count_test(P, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::range_count_test(P, querys, cnt);
    }
    else if (task == "range-report") {
        string report_qry_file = cmd.getOptionValue("-r", "range_report.qry");
        auto q_tuple = geobase::read_range_query(report_qry_file, 8, mvq::Config::get().maxSize);
        auto querys = std::get<1>(q_tuple);
        auto cnt = std::get<0>(q_tuple);
        if (algo == "mvq" || algo == "combined") ZDTest::range_report_test(P, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::range_report_test(P, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "rlog" || algo == "combined") RlogTest::range_report_test(P, querys, cnt);
        if (algo == "combined") line_splitter();
        if (algo == "pkdtree" || algo == "combined") PKDTest::range_report_test(P, querys, cnt);
    }
    else if (task == "knn") {
        size_t k = cmd.getOptionIntValue("-k", 10);
        size_t q_num = cmd.getOptionIntValue("-qn", 50000);
        if (algo == "mvq" || algo == "combined") ZDTest::knn_test(P, k, q_num);
        if (algo == "combined") line_splitter();
        if (algo == "pacz" || algo == "combined") PACZ::knn_test(P, k, q_num);
    }
    else {
        cout << "[ERROR]: Unknown task: " << task << endl;
    }
}

int main(int argc, char** argv) {
    run(argc, argv);
    return 0;
}
