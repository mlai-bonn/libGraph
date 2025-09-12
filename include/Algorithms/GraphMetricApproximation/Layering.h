//
// Created by florian on 19.10.23.
//

#ifndef TESTGRAPHLIB_LAYERING_H
#define TESTGRAPHLIB_LAYERING_H

struct LayeringCluster {
    std::vector<std::vector<NodeId>> clusters;
};

/* Creates a layering partition of the _graph with the given root node id
 *
 */
struct LayeringPartition {
    LayeringPartition() = default;
    explicit LayeringPartition(const GraphStruct& graphStruct, NodeId root_node_id, int seed = 0);

    // get root node id
    [[nodiscard]] NodeId root_node_id() const {
        return partitions[0].clusters[0][0];
    }
    std::vector<INDEX> distances;
    std::vector<LayeringCluster> partitions;
};

LayeringPartition::LayeringPartition(const GraphStruct &graph, NodeId root_node_id, int seed) {
    // create the layering tree
    std::mt19937_64 generator(seed);
    if (root_node_id == -1) {
        root_node_id = std::uniform_int_distribution<NodeId>(0, graph.nodes() - 1)(generator);
    }
    GraphStruct::BFSDistances(graph, root_node_id, distances);

    std::vector<std::vector<NodeId>> spheres;
    std::vector<GraphStruct> sphereGraphs;
    std::vector<NodeId> maxRealNodeId;
    std::vector<std::vector<NodeId>> sphereNodeIdToNodeId;
    std::vector<NodeId> nodeIdToSphereNodeId = std::vector<NodeId>(graph.nodes(), 0);
    for (INDEX i = 0; i < distances.size(); ++i) {
        // get the distance
        INDEX distance = distances[i];
        while (spheres.size() <= distance) {
            spheres.emplace_back();
            this->partitions.emplace_back();
            maxRealNodeId.emplace_back(-1);
            sphereNodeIdToNodeId.emplace_back();
        }
        spheres[distance].push_back(i);
        sphereNodeIdToNodeId[distance].push_back(i);
        nodeIdToSphereNodeId[i] = spheres[distance].size() - 1;
        ++maxRealNodeId[distance];
    }



    // create the sphere graphs
    for (int level = 0; level < spheres.size(); ++level) {
        // create the sphere _graph
        sphereGraphs.emplace_back(GraphStruct(spheres[level].size(), {}));
        for (NodeId i = 0; i < spheres[level].size(); ++i) {
            NodeId currentNodeId = spheres[level][i];
            // iterate over the neighbors of the current node
            for (NodeId j = 0; j < graph.degree(currentNodeId); ++j) {
                NodeId neighbor = graph.neighbor(currentNodeId, j);
                // check if the neighbor is in the sphere
                if (distances[neighbor] == level && i < nodeIdToSphereNodeId[neighbor]) {
                    sphereGraphs.back().AddEdge(i, nodeIdToSphereNodeId[neighbor]);
                }
            }
        }
    }


    // get the connected components of the sphere graphs
    for (int level = (int) spheres.size() - 2; level >= 0; --level) {
        auto & sphereGraph = sphereGraphs[level + 1];
        // get connected components of the sphere
        std::vector<std::vector<NodeId>> connectedComponents;
        // bfs on the sphere
        std::vector<bool> visited(sphereGraph.nodes(), false);
        for (int startNode = 0; startNode < sphereGraph.nodes(); ++startNode) {
            if (!visited[startNode]) {
                visited[startNode] = true;
                // define a queue for the current backward bfs search
                std::deque<NodeId> deque = std::deque<NodeId>();
                deque.push_back(startNode);
                // count components
                connectedComponents.emplace_back();
                connectedComponents.back().emplace_back(startNode);
                // run the backward bfs search
                while (!deque.empty()) {
                    NodeId currentNode = deque.back();
                    deque.pop_back();
                    for (NodeId i = 0; i < sphereGraph.degree(currentNode); ++i) {
                        NodeId neighbor = sphereGraph.neighbor(currentNode, i);
                        if (!visited[neighbor]) {
                            visited[neighbor] = true;
                            INDEX neighborDistance = 0;
                            if (neighbor > maxRealNodeId[level + 1]) {
                                neighborDistance = level + 1;
                            } else {
                                neighborDistance = distances[sphereNodeIdToNodeId[level + 1][neighbor]];
                            }
                            // ignore connected components
                            if (neighborDistance >= level + 1) {
                                deque.push_front(neighbor);
                                connectedComponents.back().emplace_back(neighbor);
                            }
                        }
                    }
                }
            }
        }
        // add the connected components to the partition
        // iterate over the connected components
        for (auto & connectedComponent : connectedComponents) {
            this->partitions[level + 1].clusters.emplace_back();
            for (auto & node : connectedComponent) {
                if (node <= maxRealNodeId[level + 1]) {
                    NodeId graphNodeId = sphereNodeIdToNodeId[level + 1][node];
                    this->partitions[level + 1].clusters.back().emplace_back(graphNodeId);
                }
            }

            GraphStruct & levelGraph = sphereGraphs[level];
            //modify sphere _graph of the current level
            NodeId id = levelGraph.add_node();

            std::vector<bool> levelEdges = std::vector<bool>(levelGraph.nodes(), false);
            for (auto & node : connectedComponent) {
                if (node <= maxRealNodeId[level + 1]) {
                    NodeId graphNodeId = sphereNodeIdToNodeId[level + 1][node];
                    for (NodeId i = 0; i < graph.degree(graphNodeId); ++i) {
                        NodeId neighbor = graph.neighbor(graphNodeId, i);
                        if (distances[neighbor] == level && !levelEdges[nodeIdToSphereNodeId[neighbor]]) {
                            levelGraph.AddEdge(id, nodeIdToSphereNodeId[neighbor], false);
                            levelEdges[nodeIdToSphereNodeId[neighbor]] = true;
                        }
                    }
                }
            }
        }
    }
    // add the root node to the first partition
    this->partitions[0].clusters.emplace_back();
    this->partitions[0].clusters.back().emplace_back(root_node_id);
}

/* Creates a layering tree of the _graph with the given root node id based on a _seed value if the root node is not given a random one is chosen
 *
 */

struct LayeringTree {
    LayeringTree() = default;
    explicit LayeringTree(GraphStruct& graph, int seed = 0, NodeId root_node_id = -1);
    LayeringPartition layeringPartition;
    GraphStruct tree;
};

LayeringTree::LayeringTree(GraphStruct &graph, int seed, NodeId root_node_id) : layeringPartition(graph, root_node_id, seed) {
    // create the layering tree
    std::mt19937_64 generator(seed);
    tree = GraphStruct(graph.nodes(), {});
    // iterate over the layering partitions
    for (int i = 1; i < layeringPartition.partitions.size(); ++i) {
        // iterate over the clusters of the current layering partition
        for (auto &cluster: layeringPartition.partitions[i].clusters) {
            int rand_min_pos = 0;
            // get random element from the cluster
            int rand_pos = std::uniform_int_distribution<int>(rand_min_pos, ((int) cluster.size()) - 1)(generator);
            NodeId rand_node = cluster[rand_pos];
            std::swap(cluster[rand_pos], cluster[rand_min_pos]);
            ++rand_min_pos;
            // check if the rand_node has a neighbor in the previous layer
            NodeId neighbor;
            for (NodeId j = 0; j < graph.degree(rand_node); ++j) {
                neighbor = graph.neighbor(rand_node, j);
                if (layeringPartition.distances[neighbor] == layeringPartition.distances[rand_node] - 1) {
                    break;
                }
            }
            // add edges from all nodes in the cluster to the neighbor
            for (auto &node : cluster) {
                tree.AddEdge(node, neighbor);
            }
        }
    }
    this->tree.SetType(GraphType::TREE);
}

#endif //TESTGRAPHLIB_LAYERING_H
