import os
import re

def process_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    original = content

    # --- Phase 1: Remove cpamz ---
    # Remove includes
    content = re.sub(r'#include\s*[<"](?:silva/index/|tests/)?cpamz\.hpp[>"]\n', '', content)
    content = re.sub(r'#include\s*[<"](?:silva/index/|tests/)?test_cpamz\.hpp[>"]\n', '', content)
    
    # In main.cpp, remove test execution
    content = re.sub(r'test_cpamz\(\);\n?', '', content)

    # In verify_bench.cpp, we don't have cpamz instantiated based on earlier grep, but just in case,
    # if there's any stray CPAMZ:: or cpamz_master, we should be careful. 
    # Our previous grep showed no `run_algo == "CPAMZ"` in verify_bench.cpp currently active.
    
    # --- Phase 2: Rename cpambb -> pacz ---
    content = content.replace('cpambb.hpp', 'pacz.hpp')
    content = content.replace('test_cpambb.hpp', 'test_pacz.hpp')
    content = content.replace('CPAMBB', 'PACZ')
    content = content.replace('cpambb', 'pacz')
    
    # --- Phase 3: Rename mvzd -> mvq ---
    content = content.replace('MVZD', 'MVQ')
    content = content.replace('mvzd', 'mvq')

    if content != original:
        with open(filepath, 'w') as f:
            f.write(content)
        print(f"Updated: {filepath}")

# Scan silva codebase
for root_dir in ['src', 'include/silva', 'benchmarks', 'tests']:
    if not os.path.exists(root_dir): continue
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith(('.cpp', '.hpp', '.h', '.c')):
                process_file(os.path.join(subdir, file))

