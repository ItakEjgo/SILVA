import sys

with open('benchmarks/verify_bench.cpp', 'r') as f:
    lines = f.readlines()

new_lines = []
skip = False
for i, line in enumerate(lines):
    if "namespace bg = boost::geometry;" in line:
        skip = True
        new_lines.append("#include <silva/baselines/rlog_tree.hpp>\n")
        new_lines.append("#include <silva/core/benchmark_utils.hpp>\n\n")
    
    if skip and "struct YearData {" in line:
        skip = False
        
    if not skip:
        new_lines.append(line)

with open('benchmarks/verify_bench.cpp', 'w') as f:
    f.writelines(new_lines)
