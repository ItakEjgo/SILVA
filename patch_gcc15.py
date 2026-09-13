import os
import re

# 1. Add #include <span>
for f in ["include/silva/index/mvq.hpp", "tests/test_mvq.hpp"]:
    with open(f, 'r') as file:
        content = file.read()
    if '#include <span>' not in content:
        with open(f, 'w') as file:
            file.write('#include <span>\n' + content)

# 2. Add -fpermissive to CMakeLists.txt
f = "CMakeLists.txt"
with open(f, 'r') as file:
    content = file.read()
if '-fpermissive' not in content:
    content = content.replace('set(CMAKE_CXX_STANDARD 20)', 'set(CMAKE_CXX_STANDARD 20)\nset(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")')
    with open(f, 'w') as file:
        file.write(content)

# 3. Fix time_loop.h inline issue
f = "third_party/helper/time_loop.h"
with open(f, 'r') as file:
    content = file.read()
content = content.replace('double\ntime_loop', 'inline double\ntime_loop')
with open(f, 'w') as file:
    file.write(content)

# 4. Fix std::span conversions dynamically
f = "tests/test_mvq.hpp"
with open(f, 'r') as file:
    content = file.read()

content = re.sub(r'zdtree\.build\(([^)]+)\)', r'zdtree.build(std::span<geobase::Point>(\1.data(), \1.size()))', content)
content = re.sub(r'zdtree\.multi_version_batch_insert_sorted\(([^,]+),\s*([^)]+)\)', r'zdtree.multi_version_batch_insert_sorted(std::span<geobase::Point>(\1.data(), \1.size()), \2)', content)
content = re.sub(r'zdtree\.multi_version_batch_delete_sorted\(([^,]+),\s*([^)]+)\)', r'zdtree.multi_version_batch_delete_sorted(std::span<geobase::Point>(\1.data(), \1.size()), \2)', content)

# Clean up any double-applications from previous sed commands
content = content.replace('.data().data()', '.data()').replace('.size().size()', '.size()')
content = content.replace('std::span<geobase::Point>(std::span<geobase::Point>(', 'std::span<geobase::Point>(')

with open(f, 'w') as file:
    file.write(content)

print("GCC 15 Python Patch Applied Successfully!")
