#!/bin/bash

echo "Applying patches for GCC 15 strictness..."

# 1. 修复 PKD Tree 子模块的 typename 缺失问题
sed -i 's/inline ParallelKDtree<point>::dim_type/inline typename ParallelKDtree<point>::dim_type/g' baselines/pkdtree/include/cpdd/utility/dimensinality.hpp
sed -i 's/using coord = point::coord/using coord = typename point::coord/g' baselines/pkdtree/include/cpdd/utility/tree_node.hpp
sed -i 's/using coords = point::coords/using coords = typename point::coords/g' baselines/pkdtree/include/cpdd/utility/tree_node.hpp
sed -i 's/using coord = point::coord/using coord = typename point::coord/g' baselines/pkdtree/include/cpdd/query_op/nn_search_helpers.h

# 2. 补充子模块遗漏的 std::ranges 和 algorithm 头文件
sed -i '1i #include <ranges>\n#include <algorithm>' baselines/pkdtree/include/cpdd/batch_op/build_tree.hpp
sed -i '1i #include <ranges>\n#include <algorithm>' baselines/pkdtree/include/cpdd/batch_op/batch_delete.hpp

# 3. 修复 GCC 15 对 std::span 的隐式转换严苛审查
sed -i -E 's/(zdtree\.build\()([A-Za-z0-9_]+)(\))/\1std::span<geobase::Point>(\2.data(), \2.size())\3/g' tests/test_mvq.hpp
sed -i -E 's/(zdtree\.multi_version_batch_insert_sorted\()([A-Za-z0-9_]+)(,)/\1std::span<geobase::Point>(\2.data(), \2.size())\3/g' tests/test_mvq.hpp
sed -i -E 's/(zdtree\.multi_version_batch_delete_sorted\()([A-Za-z0-9_]+)(,)/\1std::span<geobase::Point>(\2.data(), \2.size())\3/g' tests/test_mvq.hpp

echo "Patching complete! You can now run make -j in the build directory."
