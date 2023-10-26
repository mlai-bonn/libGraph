# _graph-lib: A c++ library for different _graph based tasks

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

## Graph classes

GraphStruct: undirected (labeled) _graph
DGraphStruct: directed (labeled) _graph
DDGraphStruct: directed _graph with node and edge features

GraphData<T>: collection of graphs from the above classes

## Algorithms

### Hops
1. Construct class: ```Hops(GraphData, Parameters)```
2. Run hops: ```Hops::Run(int graphId, GraphStruct pattern, RunParameters)```





