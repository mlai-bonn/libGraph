//
// Created by florian on 23.10.23.
//

#ifndef TESTGRAPHLIB_SPANNINGTREES_H
#define TESTGRAPHLIB_SPANNINGTREES_H


#include <stack>

/**
 * @brief Generates a breadth first search spanning tree of the given _graph.
 * @param graph The _graph to generate the spanning tree from.
 * @param tree The resulting spanning tree.
 * @param root_node_id The root node of the spanning tree.
 * @param deterministic If true, the algorithm will use the neighbors in the order they are stored in the _graph, otherwise it will use a random order. The random order is a little bit slower. If the order does not play a role deterministic = true should be used.
 * @param seed The seed for the random number generator. Only used if deterministic = false.
 * @return The number of components in the _graph.
 */
static int BFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, std::vector<bool>& visited, std::vector<INDEX>& distances, bool deterministic = true, int seed = 0);

/**
 * @brief Generates a depth first search spanning tree of the given _graph.
 * @param graph The _graph to generate the spanning tree from.
 * @param tree The resulting spanning tree.
 * @param root_node_id The root node of the spanning tree.
 * @param deterministic If true, the algorithm will use the neighbors in the order they are stored in the _graph, otherwise it will use a random order. The random order is a little bit slower. If the order does not play a role deterministic = true should be used.
 * @param seed The seed for the random number generator. Only used if deterministic = false.
  * @return The number of components in the _graph.
 */
static int DFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, std::vector<bool>& visited, std::vector<INDEX>& distances, bool deterministic = true, int seed = 0);



inline int BFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, std::vector<bool>& visited, std::vector<INDEX>& distances, bool deterministic, int seed){
    int components = 0;
    tree = GraphStruct(graph.nodes(), {});
    if (graph.nodes() > 0) {
        ++components;


        if (visited.size() != graph.nodes()) {
            visited.resize(graph.nodes(), false);
            std::fill(visited.begin(), visited.end(), false);
        } else {
            std::fill(visited.begin(), visited.end(), false);
        }
        if (distances.size() != graph.nodes()) {
            distances.resize(graph.nodes(), -1);
            std::fill(distances.begin(), distances.end(), -1);
        } else {
            std::fill(distances.begin(), distances.end(), -1);
        }


        std::mt19937_64 generator(seed);
        tree = GraphStruct(graph.nodes(), {});
        std::queue<NodeId> queue;
        INDEX number_nodes_visited = 0;
        std::vector<NodeId> unvisited_nodes;
        queue.push(root_node_id);
        distances[root_node_id] = 0;
        visited[root_node_id] = true;


        std::vector<NodeId> non_deterministic = std::vector<NodeId>();
        std::deque<std::pair<NodeId, NodeId>> swap_pairs;
        if (!deterministic) {
            non_deterministic.resize(graph.maxDegree, 0);
            std::iota(non_deterministic.begin(), non_deterministic.end(), 0);
        }

        while (!queue.empty() || number_nodes_visited < graph.nodes()) {
            if (queue.empty()) {
                if (components == 1) {
                    for (int i = 0; i < graph.nodes(); ++i) {
                        if (!visited[i]) {
                            unvisited_nodes.push_back(i);
                        }
                    }
                }
                ++components;
                while (visited[unvisited_nodes.back()]) {
                    unvisited_nodes.pop_back();
                }
                if (!unvisited_nodes.empty()) {
                    queue.push(unvisited_nodes.back());
                    distances[unvisited_nodes.back()] = 0;
                    visited[unvisited_nodes.back()] = true;
                } else {
                    break;
                }
            } else {
                NodeId current_node = queue.front();
                if (visited[current_node]) {
                    ++number_nodes_visited;
                }
                queue.pop();
                if (deterministic) {
                    for (auto neighbor: graph.get_neighbors(current_node)) {
                        if (!visited[neighbor]) {
                            visited[neighbor] = true;
                            queue.push(neighbor);
                            distances[neighbor] = distances[current_node] + 1;
                            tree.add_edge_no_check(current_node, neighbor);
                        }
                    }
                } else {
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
                            distances[neighbor] = distances[current_node] + 1;
                            tree.add_edge_no_check(current_node, neighbor);
                        }
                    }
                    // undo the swaps
                    while (!swap_pairs.empty()) {
                        std::pair<NodeId, NodeId> swap_pair = swap_pairs.back();
                        std::swap(non_deterministic[swap_pair.first], non_deterministic[swap_pair.second]);
                        swap_pairs.pop_back();
                    }
                }
            }
        }
    }
    if (components == 1) {
        tree.set_type(GraphType::TREE);
    }
    return components;
}

inline int DFSSpanningTree(const GraphStruct& graph, GraphStruct& tree, NodeId root_node_id, std::vector<bool>& visited, std::vector<INDEX>& distances, bool deterministic, int seed){
    int components = 0;
    tree = GraphStruct(graph.nodes(), {});
    if (graph.nodes() > 0) {
        ++components;

        if (visited.size() != graph.nodes()) {
            visited.resize(graph.nodes(), false);
            std::fill(visited.begin(), visited.end(), false);
        } else {
            std::fill(visited.begin(), visited.end(), false);
        }
        if (distances.size() != graph.nodes()) {
            distances.resize(graph.nodes(), -1);
            std::fill(distances.begin(), distances.end(), -1);
        } else {
            std::fill(distances.begin(), distances.end(), -1);
        }

        std::mt19937_64 generator(seed);
        tree = GraphStruct(graph.nodes(), {});
        std::stack<NodeId> stack;
        INDEX number_nodes_visited = 0;
        std::vector<NodeId> unvisited_nodes;
        stack.push(root_node_id);
        distances[root_node_id] = 0;
        visited[root_node_id] = true;

        std::vector<NodeId> non_deterministic = std::vector<NodeId>();
        std::deque<std::pair<NodeId, NodeId>> swap_pairs;
        if (!deterministic) {
            non_deterministic.resize(graph.maxDegree, 0);
            std::iota(non_deterministic.begin(), non_deterministic.end(), 0);
        }

        while (!stack.empty() || number_nodes_visited < graph.nodes()) {
            if (stack.empty()) {
                if (components == 1) {
                    for (int i = 0; i < graph.nodes(); ++i) {
                        if (!visited[i]) {
                            unvisited_nodes.push_back(i);
                        }
                    }
                }
                ++components;
                while (visited[unvisited_nodes.back()]) {
                    unvisited_nodes.pop_back();
                }
                if (!unvisited_nodes.empty()) {
                    stack.push(unvisited_nodes.back());
                    distances[unvisited_nodes.back()] = 0;
                    visited[unvisited_nodes.back()] = true;
                } else {
                    break;
                }
            } else {
                NodeId current_node = stack.top();
                stack.pop();
                if (visited[current_node]) {
                    ++number_nodes_visited;
                }
                if (deterministic) {
                    for (auto neighbor: graph.get_neighbors(current_node)) {
                        if (!visited[neighbor]) {
                            stack.push(neighbor);
                            visited[neighbor] = true;
                            tree.add_edge_no_check(current_node, neighbor);
                            distances[neighbor] = distances[current_node] + 1;
                        }
                    }
                } else {
                    // iterate randomly over the neighbors of the current node
                    INDEX degree = graph.degree(current_node);
                    for (NodeId i = 0; i < degree; ++i) {
                        NodeId rand_idx = std::uniform_int_distribution<NodeId>(i, degree - 1)(generator);
                        NodeId neighbor = graph.get_neighbors(current_node)[non_deterministic[rand_idx]];
                        std::swap(non_deterministic[rand_idx], non_deterministic[i]);
                        swap_pairs.emplace_back(rand_idx, i);
                        if (!visited[neighbor]) {
                            stack.push(neighbor);
                            visited[neighbor] = true;
                            tree.add_edge_no_check(current_node, neighbor);
                            distances[neighbor] = distances[current_node] + 1;
                        }
                    }
                    // undo the swaps
                    while (!swap_pairs.empty()) {
                        std::pair<NodeId, NodeId> swap_pair = swap_pairs.back();
                        std::swap(non_deterministic[swap_pair.first], non_deterministic[swap_pair.second]);
                        swap_pairs.pop_back();
                    }
                }
            }
        }
    }
    if (components == 1) {
        tree.set_type(GraphType::TREE);
    }
    return components;
}


#endif //TESTGRAPHLIB_SPANNINGTREES_H
