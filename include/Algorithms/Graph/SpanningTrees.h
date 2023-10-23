//
// Created by florian on 23.10.23.
//

#ifndef TESTGRAPHLIB_SPANNINGTREES_H
#define TESTGRAPHLIB_SPANNINGTREES_H


#include <stack>

/**
 * @brief Generates a breadth first search spanning tree of the given graph.
 * @param graph The graph to generate the spanning tree from.
 * @param tree The resulting spanning tree.
 * @param root_node_id The root node of the spanning tree.
 * @param deterministic If true, the algorithm will use the neighbors in the order they are stored in the graph, otherwise it will use a random order. The random order is a little bit slower. If the order does not play a role deterministic = true should be used.
 * @param seed The seed for the random number generator. Only used if deterministic = false.
 */
static void BFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, bool deterministic = true, int seed = 0);

/**
 * @brief Generates a depth first search spanning tree of the given graph.
 * @param graph The graph to generate the spanning tree from.
 * @param tree The resulting spanning tree.
 * @param root_node_id The root node of the spanning tree.
 * @param deterministic If true, the algorithm will use the neighbors in the order they are stored in the graph, otherwise it will use a random order. The random order is a little bit slower. If the order does not play a role deterministic = true should be used.
 * @param seed The seed for the random number generator. Only used if deterministic = false.
 */
static void DFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, bool deterministic = true, int seed = 0);



inline void BFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, bool deterministic, int seed){
    std::mt19937_64 generator(seed);
    tree = GraphStruct(graph.nodes(), {});
    std::queue<NodeId> queue;
    std::vector<bool> visited(graph.nodes(), false);
    queue.push(root_node_id);

    std::vector<NodeId> non_deterministic = std::vector<NodeId>();
    std::deque<std::pair<NodeId, NodeId>> swap_pairs;
    if (!deterministic){
        non_deterministic.resize(graph.maxDegree, 0);
        std::iota(non_deterministic.begin(), non_deterministic.end(), 0);
    }

    while (!queue.empty()){
        NodeId current_node = queue.front();
        visited[current_node] = true;
        queue.pop();
        if (deterministic){
            for (auto neighbor : graph.get_neighbors(current_node)) {
                if (!visited[neighbor]) {
                    visited[current_node] = true;
                    queue.push(neighbor);
                    tree.add_edge_linear(current_node, neighbor);
                }
            }
        }
        else{
            // iterate randomly over the neighbors of the current node
            INDEX degree = graph.degree(current_node);
            for (NodeId i = 0; i < degree; ++i) {
                NodeId rand_idx = std::uniform_int_distribution<NodeId>(i, degree - 1)(generator);
                NodeId neighbor = graph.get_neighbors(current_node)[non_deterministic[rand_idx]];
                std::swap(non_deterministic[rand_idx], non_deterministic[i]);
                swap_pairs.emplace_back(rand_idx, i);
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                    tree.add_edge_linear(current_node, neighbor);
                }
            }
            // undo the swaps
            while (!swap_pairs.empty()){
                std::pair<NodeId, NodeId> swap_pair = swap_pairs.back();
                std::swap(non_deterministic[swap_pair.first], non_deterministic[swap_pair.second]);
                swap_pairs.pop_back();
            }
        }
    }
    tree.isTree = true;
    tree.set_type(GraphType::TREE);
}

inline void DFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, bool deterministic, int seed){
    std::mt19937_64 generator(seed);
    tree = GraphStruct(graph.nodes(), {});
    std::stack<NodeId> queue;
    std::vector<bool> visited(graph.nodes(), false);
    queue.push(root_node_id);

    std::vector<NodeId> non_deterministic = std::vector<NodeId>();
    std::deque<std::pair<NodeId, NodeId>> swap_pairs;
    if (!deterministic){
        non_deterministic.resize(graph.maxDegree, 0);
        std::iota(non_deterministic.begin(), non_deterministic.end(), 0);
    }

    while (!queue.empty()){
        NodeId current_node = queue.top();
        visited[current_node] = true;
        queue.pop();
        if (deterministic){
            for (auto neighbor : graph.get_neighbors(current_node)) {
                if (!visited[neighbor]) {
                    visited[current_node] = true;
                    queue.push(neighbor);
                    tree.add_edge_linear(current_node, neighbor);
                }
            }
        }
        else{
            // iterate randomly over the neighbors of the current node
            INDEX degree = graph.degree(current_node);
            for (NodeId i = 0; i < degree; ++i) {
                NodeId rand_idx = std::uniform_int_distribution<NodeId>(i, degree - 1)(generator);
                NodeId neighbor = graph.get_neighbors(current_node)[non_deterministic[rand_idx]];
                std::swap(non_deterministic[rand_idx], non_deterministic[i]);
                swap_pairs.emplace_back(rand_idx, i);
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                    tree.add_edge_linear(current_node, neighbor);
                }
            }
            // undo the swaps
            while (!swap_pairs.empty()){
                std::pair<NodeId, NodeId> swap_pair = swap_pairs.back();
                std::swap(non_deterministic[swap_pair.first], non_deterministic[swap_pair.second]);
                swap_pairs.pop_back();
            }
        }
    }
    tree.isTree = true;
    tree.set_type(GraphType::TREE);
}


#endif //TESTGRAPHLIB_SPANNINGTREES_H
