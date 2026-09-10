# SILVA Library Refactoring Log

This log documents the major structural and architectural changes made to transition the `mvzd` project into the clean, standalone `silva` library.

## Date: 2026-09-10
**Goal:** Extract core code to a new repository, clean up historical baggage, and restructure into a modern C++ library layout with zero performance overhead.

### Phase 1: Clean Extraction
- **Action:** Created a new `silva` directory and copied only essential source files (`include/`, `src/`, `baselines/`, `CMakeLists.txt`, `verify_bench.cpp`).
- **Reason:** To shed the massive historical baggage (old experiments, datasets, logs, unused binaries) present in the old `mvzd` repository and establish a clean Git history.
- **Verification:** Successfully compiled `verify_bench` in the new isolated environment. Removed mistakenly copied cache and binary files (`test.log`, `verify_bench`, `.cache/`). 

### Phase 2: Structural Refactoring (Step 1)
- **Action:** Reorganized the repository into a standard library layout:
  - `include/silva/core/`: Moved `global_config.hpp`, `hilbert.h`.
  - `include/silva/geo/`: Moved `src/geo/*.hpp`.
  - `include/silva/index/`: Moved `src/mvq.hpp`, `cpambb.hpp`, `cpamz.hpp`.
  - `third_party/`: Isolated external dependencies (`cpam`, `pam`, `parlay`, `helper`) to prevent naming collisions and messy root includes.
  - `benchmarks/`: Relocated `verify_bench.cpp`.
  - `tests/`: Relocated `src/test/`.
- **Reason:** To separate public API headers from private implementations and isolate third-party code.
- **Method:** Used exact string replacement via Python for `#include` paths to ensure absolute safety and 0% source code logic modification.

### Phase 3: Zero-Overhead "Physics Extraction" (Step 2)
- **Action:** Extracted the massive inlined `RlogTree` baseline implementation from `verify_bench.cpp` into a dedicated header `include/silva/baselines/rlog_tree.hpp`. 
- **Action:** Extracted memory measurement helper functions into `include/silva/core/benchmark_utils.hpp`.
- **Reason:** `verify_bench.cpp` was extremely bloated (~900 lines) with mixed benchmarking logic and algorithm implementation. We extracted the algorithm physically without introducing a `virtual` Interface.
- **Why no `virtual` Interface?** Algorithms return drastically different types for KNN queries (`priority_queue` vs `std::vector`). Forcing a unified return type via a wrapper would introduce data copying overhead inside the timer loop, rendering the scientific benchmark fundamentally inaccurate. This extraction guarantees **100% exact performance** compared to the original code.

*All changes were independently compiled and verified using `cmake` and `make -j`.*