import sys

def comment_out(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    for line in lines:
        if 'CPAMZ::' in line or 'CPAMZ_init' in line or 'zMAP_init' in line or 'zMAP_range_report' in line:
            new_lines.append('// ' + line)
        else:
            new_lines.append(line)
            
    with open(filepath, 'w') as f:
        f.writelines(new_lines)

comment_out('src/main.cpp')
comment_out('tests/test_mvq.hpp')
