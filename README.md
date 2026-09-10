# SILVA: Spatial Index Library with Version Access

## Overview
**SILVA** is a high-performance main-memory spatial index library for managing multi-version spatial data (points and polygons). It enables storing multiple versions of spatial indexes and supports real-time spatial queries on any stored version, including range counting/reporting and k-nearest neighbors (KNN). 
In addition to *standard spatial queries*, it also supports
*version control* operations (commit, merge, spatial diff, purge) are optimized for parallel execution. 
Our experimental results show that SILVA can achieve  **100x faster operations** and **70% lower memory usage** than state-of-the-art multi-version spatial index MV(3)R-tree from *libspatialindex* library.

**Ploygon Support** SILVA also supports processing polygons like traditional *R-trees*. We also integrate rectangle intersects queries, which is widely used in real-world applications.

**Disk Interaction** This capability is particularly helpful when certain versions need to be archived for long-term importance, or when accumulated versions cause significant high memory pressure. In such scenarios, specific versions can be stored from memory to disk while keeping the calculated information. This ensures when the versions are later retrieved from disk, no extra computation overhead is required to reconstruct the index.

## Key Features
- **Multi-Version Control**  
   versioned operations:  
   - `commit`
   - `merge`
   - `spatial diff`, and 
   - `purge`
- **Standrad Single Version Spatial Queries**  
  Execute on any historical version:
  - `range_count` / `range_report`
  - `knn_search`
- **Polygon Support**
  - Build, update, and intersect query for rectangles.
- **Disk Interactions**
  - Store/Load an entire version to/from disk.
- **Parallel Updates**
  - Update to the multi-version index can be highly in parallel.  
- **High Efficiency**  
  - Single-version queries match traditional spatial indexes (Rtrees in boost) with no extra overhead.

---

## Implemented Indexes
### 1. PaC-Z tree
- **High-level Idea**  
  Maps spatial points to 1D via Z-curve, then extends the **PaC Tree** (1D functional binary search tree) for versioned operations and spatial queries.

### 2. MVQ Tree (Multi-version Quadtree)
- **High-level Idea**  
  Leverages **spatial invariance** to maximize sharing of unchanged regions across versions.
  Achieves minimal memory footprint as well as fast versioned operations.
---

## Performance Comparison
- Please refer to Experiment section in our paper and supplementary material.
---

## Code Structure
./

├── src/ # SILVA core

├── include/ # PaC-Tree libraries

Other files.
The implementation of polygon support and disk interaction can be found in branch *silva-box*.

## Compilation
bash

cd src/ # or baselines/

make # Generates binaries

---

## Command-Line Parameters
| Parameter              | Description                          | Example                   |
|------------------------|--------------------------------------|---------------------------|
| `-i <Path-to-Input>`   | Input data file path                 | `-i data/points.in`      |
| `-o <Path-to-Output>`  | Output results path                  | `-o results.log`      |
| `-t <Task-Name>`       | Operation type (`build`, `insert`,`range-count`, `box-intersect` etc.)| `-t build`          |
| `-a <Algorithm-Name>` | Index type (`MVQ`, `paczz`, `paczbb`, etc.) | `-a MVQ`                |
| `-b <Batch-file>`      | Batch operations file                | `-b batch/ops.in`        |
| `-bf <batch-fraction>` | Percent of batch to process | `-bf 10`               |
| `-real <0/1>`   | Real-world dataset flag              | `-real 1`              |

### Example Command
bash

./main -i data/uniform_100M.in -o build.log -t build -a MVQ -real 0
