#include <iostream>
#include <string>
#include <cpam/parse_command_line.h>

#include <silva/core/dataset_loader.hpp>
#include "benchmark_engine.hpp"
#include "runners/runner_common.hpp"
#include "runners/mvq_runner.hpp"
#include "runners/pacz_runner.hpp"
#include "runners/rlog_runner.hpp"
#include "runners/pkd_log_runner.hpp"

int main(int argc, char** argv) {
    cpam::commandLine cmd(argc, argv, "");
    std::string run_algo = cmd.getOptionValue("-algo", "all");
    std::string dataset_dir = cmd.getOptionValue("-dir", "/data/bhuan102/SILVA-dataset/bhutan_workload");
    int start_year = cmd.getOptionIntValue("-start_year", 2018);
    int end_year = cmd.getOptionIntValue("-end_year", 2026);
    size_t q_step = cmd.getOptionIntValue("-q_step", 1000);
    std::string commits_dir = cmd.getOptionValue("-commits_dir", "01_commits");

    // 1. In-Memory Data Preparation
    DatasetLoader data;
    data.load(dataset_dir, commits_dir, start_year, end_year);

    // 2. Base Data Conversion for Runners
    std::vector<Value> P_base;
    P_base.reserve(data.base_data.size());
    for(const auto& pt : data.base_data) {
        P_base.push_back(std::make_pair(BoostPoint(pt.x, pt.y), pt.id));
    }

    // 3. Execution Engine
    if (run_algo == "MVQ" || run_algo == "all") {
        MVQRunner runner(mvq::Config::get().leaf_size);
        runner.build_base(P_base);
        BenchmarkEngine<MVQRunner>::run_workload(runner, data, "MVQ", q_step);
    }
    if (run_algo == "PACZ" || run_algo == "all") {
        PACZRunner runner;
        runner.build_base(P_base);
        BenchmarkEngine<PACZRunner>::run_workload(runner, data, "PACZ", q_step);
    }
    if (run_algo.find("Rlog") != std::string::npos || run_algo == "all") {
        // Simplified to run 1-year compaction as an example
        RlogRunner runner(1);
        runner.build_base(P_base);
        BenchmarkEngine<RlogRunner>::run_workload(runner, data, "Rlog_1yr", q_step);
    }

        if (run_algo.find("Pkd") != std::string::npos || run_algo == "all" || run_algo == "combined") {
        PKDLogRunner runner;
        runner.build_base(P_base);
        BenchmarkEngine<PKDLogRunner>::run_workload(runner, data, "PkdLogTree", q_step);
    }

    return 0;
}
