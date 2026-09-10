#pragma once
#include <atomic>
#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <silva/index/pacz.hpp>
#include <silva/index/mvq.hpp>
#include <silva/baselines/rlog_tree.hpp>

struct BranchRes { double fork_ms; double commit_ms; double merge_ms; double mem_mb; };

inline std::pair<double, double> mem_pacz(const PACZ::zmap& latest_branch) {
    return {cpam::cpam_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0), 0.0};
}

inline std::pair<double, double> mem_mvq() { 
    return {mvq::global_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0), 0.0}; 
}

inline std::pair<double, double> mem_rlog(const RlogTree& master, const std::vector<RlogBranch>& history) {
    double rtree_mem = boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
    double delta_mem = 0;
    for (const auto& b : history) {
        delta_mem += (b.insert_log.capacity() + b.remove_log.capacity()) * sizeof(Value);
    }
    delta_mem += (master.insert_log.capacity() + master.remove_log.capacity()) * sizeof(Value);
    delta_mem += unordered_set_mem_estimate(master.removed_ids.bucket_count()) + master.removed_ids.size() * sizeof(size_t);
    return {rtree_mem, delta_mem / (1024.0 * 1024.0)};
}

inline void log_branch(std::ofstream& fout, const std::string& algo, int batch_idx, BranchRes res) {
    fout << std::flush;
    fout << algo << " | " << batch_idx << " | " << res.fork_ms << " | " << res.commit_ms << " | " << res.merge_ms << " | " << res.mem_mb << "\n";
}
