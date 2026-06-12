# libGraph: A C++ library for graph-based tasks

libGraph is a header-only C++20 library providing graph data structures, binary graph serialization, and implementations/wrappers of several graph algorithms (HOPS pattern counting, geodesic cores, graph Tukey depth, outerplanar subgraphs, spanning trees, and graph edit distance via [gedlib](https://github.com/dbblumenthal/gedlib)).

Everything is included via the umbrella header [include/libGraph.h](include/libGraph.h).

## Repository layout

| Path | Content |
|---|---|
| `include/` | The library itself (header-only) |
| `include/GraphDataStructures/` | Graph classes (`GraphBase.h`, `GraphData.h`, …) |
| `include/Algorithms/` | GED wrapper, geodesic core, outerplanar subgraphs, spanning trees, metric approximation |
| `include/Hops/` | HOPS embedding/pattern counting |
| `include/Closures/` | Graph closure operators |
| `include/io/` | File evaluation and I/O helpers |
| `examples/` | Standalone example projects (`graphs`, `tu_datasets`, `algorithms/ged`), each with its own `CMakeLists.txt` |
| `executables/` | Experiment/utility binaries (`graph_conversions`, `test`, `geodesic_core`) |
| `Google_tests/` | GoogleTest test suite (`Google_tests/CMakeLists.txt`) |
| `python/graph_plotting/` | Python helper for plotting graphs |
| `cmake/` | `FindGUROBI.cmake` for the optional Gurobi-based exact GED methods |

Note: there is no root `CMakeLists.txt`; each example/executable/test project is configured separately and adds `include/` to its include path.

## Requirements

- C++20 compiler (some example projects pin `g++-11`)
- CMake ≥ 3.20
- OpenMP (used by HOPS and the test/example projects)
- Boost Graph (used by the GED evaluation code)
- Optional: [gedlib](https://github.com/dbblumenthal/gedlib) for graph edit distance (enabled with `#define GEDLIB`)
- Optional: Gurobi for the exact GED methods F1/F2/COMPACT_MIP (enabled with `#define GUROBI`)

## Graph classes

- **GraphStruct**: undirected (labeled) graph
- **DGraphStruct**: directed (labeled) graph
- **DDGraphStruct**: directed graph with node and edge features
- **UDGraphStruct**: undirected graph with node and edge features
- **GraphData\<T\>**: collection of graphs of one of the above classes, for saving and loading multiple graphs in one file

## Algorithms

### HOPS [(paper)](https://dl.acm.org/doi/10.1145/3394486.3403180)
Estimating pattern (tree/graph) embedding counts by sampling.
1. Construct class: `Hops(GraphData, Parameters)`
2. Run hops: `Hops::Run(int graphId, GraphStruct pattern, RunParameters)`

### Geodesic Core [(paper)](https://arxiv.org/abs/2206.07350)
See `include/Algorithms/GeodesicCore/`.

### Graph Tukey Depth [(paper)](https://ceur-ws.org/Vol-3341/KDML-LWDA_2022_CRC_705.pdf)

### Outerplanar subgraphs & spanning trees
Maximal outerplanar subgraph algorithms (Mitchell, DFS-based) and spanning tree utilities in `include/Algorithms/Graph/`.

### Graph Edit Distance (GED)
For GED computations we wrap the external library [gedlib](https://github.com/dbblumenthal/gedlib) (`include/Algorithms/GED/`). On top of the raw distance, the module can reconstruct *edit paths* (sequences of intermediate graphs from a source to a target graph) with configurable operation-ordering strategies.

Build steps:

1. Navigate to the [external](external) folder and clone the gedlib repository:
```bash
git clone https://github.com/dbblumenthal/gedlib.git
```
2. Follow the [instructions](https://github.com/dbblumenthal/gedlib/blob/master/README.md) in the gedlib repository to build the library.
   Hints:
    - You need: CMake, Doxygen and OpenMP
    - If using the C++20 standard (as in this repo) replace boost 1.69.0 by boost 1.89.0 [(download)](https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.gz).
    - If everything is installed properly, use:
      ```bash
      python install.py
      ```
      to build the library.
3. Define `GEDLIB` before including `libGraph.h` (and `GUROBI` if the MIP-based exact methods are wanted).
4. [examples/algorithms/ged](examples/algorithms/ged) contains a complete example using gedlib with our graph classes.

## Tests

The GoogleTest suite lives in `Google_tests/` (graph classes, algorithms, closures, cores, outerplanar subgraphs):

```bash
cd Google_tests && cmake -B build && cmake --build build && ./build/Google_Tests_run
```

The GED module is currently not covered by the test suite (see [TODO.md](TODO.md)).

## Binary graph format (`.bgfs` for graphs with fewer than 2^32 edges, `.bgf` otherwise)

|  Parameter Type  |             Description              |
|:----------------:|:------------------------------------:|
|     **int**      |    *compatibility_format_version*    |
|     **int**      |            *graph_number*            |
|                  |      **Repeat for each graph**       |
| **unsigned int** |    string length of *graph_name*     |
|    **string**    |             *graph_name*             |
|     **enum**     |             *graph_type*             |
|    **size_t**    |            *node_number*             |
| **unsigned int** |           *node_features*            |
|                  |                                      |
|                  |  **Repeat for every node feature**   |
| **unsigned int** | string length of *node_feature_name* |
|    **string**    |         *node_feature_name*          |
|                  |                                      |
|    **size_t**    |            *edge_number*             |
| **unsigned int** |           *edge_features*            |
|                  |                                      |
|                  |  **Repeat for every edge feature**   |
| **unsigned int** | string length of *edge_feature_name* |
|    **string**    |         *edge_feature_name*          |
|                  |                                      |
|                  |      **Repeat for each graph**       |
|                  |                                      |
|                  |      **Repeat for every node**       |
|    **double**    |           *node_feature_i*           |
|                  |                                      |
|    **size_t**    |             *edge_head*              |
|    **size_t**    |             *edge_tail*              |
|                  |      **Repeat for every edge**       |
|    **double**    |           *edge_feature_i*           |

Note: the format stores raw in-memory representations (`size_t`, enum values), so files are not portable across platforms with different type sizes or endianness.

## Known issues

A code review (June 2026) identified several bugs and structural improvements; see [TODO.md](TODO.md) for the prioritized list.
