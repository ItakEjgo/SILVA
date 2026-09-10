import os

replacements = {
    # Replace external include paths
    '"cpam/': '<cpam/',
    '"pam/': '<pam/',
    '"parlay/': '<parlay/',
    '"helper/': '<helper/',
    
    # Verify_bench specific
    '"src/mvq.hpp"': '<silva/index/mvq.hpp>',
    '"src/cpambb.hpp"': '<silva/index/cpambb.hpp>',
    '"src/cpamz.hpp"': '<silva/index/cpamz.hpp>',
    '"src/global_config.hpp"': '<silva/core/global_config.hpp>',
    '"src/hilbert.h"': '<silva/core/hilbert.h>',

    # Internal core includes
    '"global_config.hpp"': '<silva/core/global_config.hpp>',
    '"../global_config.hpp"': '<silva/core/global_config.hpp>',
    
    '"hilbert.h"': '<silva/core/hilbert.h>',
    '"../hilbert.h"': '<silva/core/hilbert.h>',

    '"mvq.hpp"': '<silva/index/mvq.hpp>',
    '"../mvq.hpp"': '<silva/index/mvq.hpp>',

    '"cpambb.hpp"': '<silva/index/cpambb.hpp>',
    '"../cpambb.hpp"': '<silva/index/cpambb.hpp>',

    '"cpamz.hpp"': '<silva/index/cpamz.hpp>',
    '"../cpamz.hpp"': '<silva/index/cpamz.hpp>',

    '"geo/point.hpp"': '<silva/geo/point.hpp>',
    '"../geo/point.hpp"': '<silva/geo/point.hpp>',
    '"point.hpp"': '<silva/geo/point.hpp>',

    '"geo/operations.hpp"': '<silva/geo/operations.hpp>',
    '"../geo/operations.hpp"': '<silva/geo/operations.hpp>',
    '"operations.hpp"': '<silva/geo/operations.hpp>',

    '"geo/io.hpp"': '<silva/geo/io.hpp>',
    '"../geo/io.hpp"': '<silva/geo/io.hpp>',
    '"io.hpp"': '<silva/geo/io.hpp>',

    '"test/test_mvq.hpp"': '"test_mvq.hpp"',
    '"test/test_cpambb.hpp"': '"test_cpambb.hpp"',
    '"test/test_cpamz.hpp"': '"test_cpamz.hpp"',
    '"test_utils.hpp"': '"test_utils.hpp"',
}

directories_to_scan = ['include/silva', 'src', 'benchmarks', 'tests']

for root_dir in directories_to_scan:
    if not os.path.exists(root_dir): continue
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith(('.cpp', '.hpp', '.h', '.c')):
                filepath = os.path.join(subdir, file)
                with open(filepath, 'r') as f:
                    content = f.read()
                
                new_content = content
                for old, new in replacements.items():
                    # Only replace exact include statements to be ultra-safe
                    old_include = f'#include {old}'
                    new_include = f'#include {new}'
                    new_content = new_content.replace(old_include, new_include)

                if new_content != content:
                    with open(filepath, 'w') as f:
                        f.write(new_content)
                    print(f"Updated includes in: {filepath}")

