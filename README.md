# libgraph: A c++ library for different graph based tasks

## Graph classes

We provide the following graph classes:

**GraphStruct**: undirected (labeled) _graph
**DGraphStruct**: directed (labeled) _graph
**DDGraphStruct**: directed _graph with node and edge features

**GraphData<<T>>**: collection of graphs from the above classes, for saving and loading multiple graphs in one file.

## Algorithms

### Hops [(paper)](https://dl.acm.org/doi/10.1145/3394486.3403180)
1. Construct class: ```Hops(GraphData, Parameters)```
2. Run hops: ```Hops::Run(int graphId, GraphStruct pattern, RunParameters)```

### Geodesic Core [(paper)](https://arxiv.org/abs/2206.07350)

### Graph Tukey Depth [(paper)](https://ceur-ws.org/Vol-3341/KDML-LWDA_2022_CRC_705.pdf)

### Graph Edit Distance
For GED computations, we provide a wrapper for the external library [gedlib](https://github.com/dbblumenthal/gedlib).
The following steps are necessary to build the library:

1. Navigate to the [external](external) folder and clone the gedlib repository:
```bash
git clone https://github.com/dbblumenthal/gedlib.git
```
2. Follow the instructions in the gedlib repository to build the library.
   Hints:
    - You need: CMake, Doxygen and OpenMP
    - If using C++20 standard (as in this repo) replace boost.1.69.0 by boost.1.89.0 [(download)](https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.gz).
    - If everything is installed properly, use:
   ```bash
   python install.py
   ``` 
   to build the library.
3. In [examples/algorithms/ged](examples/algorithms/ged) we provide an example for using the gedlib library with our graph classes.




## Graph Format (.bgfs for graphs with less than 2^32 edges .bgf else)

|  Parameter Type  |             Description              |
|:----------------:|:------------------------------------:|
|     **int**      |    *compatibility_format_version*    |
|     **int**      |            *graph_number*            |
|                  |      **Repeat for each _graph**       |
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
|                  |  **Repeat for every node feature**   |
| **unsigned int** | string length of *edge_feature_name* |
|    **string**    |         *edge_feature_name*          |
|                  |                                      |
|                  |      **Repeat for each _graph**       |
|                  |                                      |
|                  |      **Repeat for every node**       |
|    **double**    |           *node_feature_i*           |
|                  |                                      |
|    **size_t**    |             *edge_head*              |
|    **size_t**    |             *edge_tail*              |
|                  |      **Repeat for every edge**       |
|    **double**    |           *edge_feature_i*           |
|                  |                                      |






