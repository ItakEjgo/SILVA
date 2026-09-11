import re

with open('benchmarks/verify_bench.cpp', 'r') as f:
    content = f.read()

def print_matches(regex_str, description):
    print(f"\n=== {description} ===")
    matches = set(re.findall(regex_str, content))
    for m in matches:
        print(m.strip())

print_matches(r'.*=\s*new\s+[a-zA-Z0-9_:]+.*', "Initialization")
print_matches(r'.*(?:\.|->)build\(.*', "Build")
print_matches(r'.*(?:\.|->)range_report\(.*', "Range Query")
print_matches(r'.*(?:\.|->)knn_report\(.*', "KNN Query")
print_matches(r'.*(?:\.|->)commit\(.*', "Commit / Update")

