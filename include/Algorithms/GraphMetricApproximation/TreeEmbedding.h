//
// Created by florian on 17.10.23.
//
// This file contains the algorithm from the paper "Constant approximation algorithms for embedding
//graph metrics into trees and outerplanar graphs" by Chepoi, Dragan, Newman, Rabinovich and Vaxès.

#ifndef LIBGRAPH_TREEEMBEDDING_H
#define LIBGRAPH_TREEEMBEDDING_H

struct TreeEmbeddingInputParameters {
    NodeId root_node_id = 0;
    int seed = 0;
};

struct LayeringCluster {
    std::vector<std::vector<NodeId>> layeringCluster;
};

struct LayeringPartition {
    NodeId rootNodeId = 0;
    std::vector<LayeringCluster> layeringPartition;
};

struct LayeringTree {
    NodeId rootNodeId = 0;
    LayeringPartition layeringPartition;
    GraphStruct tree;
};

class TreeEmbeddingAlgorithm {
public:
    explicit TreeEmbeddingAlgorithm(GraphStruct &graph) : _graph(graph) {};

    void Run(TreeEmbeddingInputParameters& inputParameters){
        LayeringTree layeringTree;
        GetLayeringTree(layeringTree, inputParameters);
    }

    LayeringPartition GetLayeringPartition(const TreeEmbeddingInputParameters& inputParameters);

    void GetLayeringTree(LayeringTree& layeringTree, const TreeEmbeddingInputParameters& inputParameters){
        GetLayeringPartition(inputParameters);
    }


private:
    // Create a layering partition of the graph
private:
    GraphStruct & _graph;

};

inline LayeringPartition TreeEmbeddingAlgorithm::GetLayeringPartition(const TreeEmbeddingInputParameters& inputParameters) {
    std::vector<INDEX> distances;
    GraphStruct::BFSDistances(_graph, inputParameters.root_node_id, distances);

    std::vector<std::vector<NodeId>> spheres;
    std::vector<GraphStruct> sphereGraphs;
    std::vector<NodeId> maxRealNode;
    std::vector<std::vector<NodeId>> sphereNodeIdToNodeId;
    std::vector<NodeId> nodeIdToSphereNodeId = std::vector<NodeId>(_graph.nodes(), 0);

    std::vector<std::vector<std::vector<NodeId>>> partition;
    std::vector<int> connectedComponentIndices = std::vector<int>(_graph.nodes(), -1);
    for (INDEX i = 0; i < distances.size(); ++i) {
        // get the distance
        INDEX distance = distances[i];
        while (spheres.size() <= distance) {
            spheres.emplace_back();
            partition.emplace_back();
            maxRealNode.emplace_back(0);
            sphereNodeIdToNodeId.emplace_back();
        }
        spheres[distance].push_back(i);
        sphereNodeIdToNodeId[distance].push_back(i);
        nodeIdToSphereNodeId[i] = spheres[distance].size() - 1;
        ++maxRealNode[distance];
    }



    // create the sphere graphs
    for (int level = 0; level < spheres.size(); ++level) {
        // create the sphere graph
        sphereGraphs.emplace_back(GraphStruct(spheres[level].size(), {}));
        for (NodeId i = 0; i < spheres[level].size(); ++i) {
            NodeId currentNodeId = spheres[level][i];
            // iterate over the neighbors of the current node
            for (NodeId j = 0; j < _graph.degree(currentNodeId); ++j) {
                NodeId neighbor = _graph.neighbor(currentNodeId, j);
                // check if the neighbor is in the sphere
                if (distances[neighbor] == level && i < nodeIdToSphereNodeId[neighbor]) {
                    sphereGraphs.back().add_edge(i, nodeIdToSphereNodeId[neighbor]);
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
                        NodeId neighbor = sphereGraph.neighbor(currentNode, sphereGraph.degree(currentNode) - i);
                        // ignore connected components
                        INDEX neighborDistance = distances[neighbor];
                        if (neighborDistance >= level) {
                            deque.push_front(neighbor);
                            if (neighborDistance == level && !visited[neighbor]){
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
            partition[level + 1].emplace_back();
            for (auto & node : connectedComponent) {
                if (node < maxRealNode[level + 1]) {
                    partition[level + 1].back().emplace_back(sphereNodeIdToNodeId[level +1][node]);
                }
                else
                {
                    break;
                }
            }

            //modify sphere graph of the current level
            NodeId id = sphereGraphs[level].add_node();


            for (auto & node : connectedComponent) {
                connectedComponentIndices[node] = (int) partition[level].size() - 1;
                for (NodeId i = 0; i < _graph.degree(sphereNodeIdToNodeId[level + 1][node]); ++i) {
                    NodeId neighbor = _graph.neighbor(sphereNodeIdToNodeId[level + 1][node], i);
                    if (distances[neighbor] == level) {
                        sphereGraph.add_edge(id, neighbor);
                    }
                }
            }
        }


    }



    std::vector<bool> visited(_graph.nodes(), false);
    // start a backward search from each of the most distant vertices of the forward search until some visited vertex is reached
    std::vector<NodeId> nodeOrder;
    for (auto & sphere : spheres) {
        for (auto & node : sphere) {
            nodeOrder.push_back(node);
        }
    }
    while (!nodeOrder.empty()) {
        auto startNode = nodeOrder.back();
        nodeOrder.pop_back();
        if (!visited[startNode]) {
            // stores the levels of the nodes and the cluster index
            std::vector<std::pair<int, std::vector<NodeId>>> inverseLocalSpheres = std::vector<std::pair<int, std::vector<NodeId>>>(distances[startNode] + 1, std::make_pair(-1, std::vector<NodeId>()));

            // define a queue for the current backward bfs search
            std::deque<NodeId> nodes = std::deque<NodeId>();
            nodes.push_back(startNode);
            // run the backward bfs search
            while (!nodes.empty()) {
                NodeId currentNode = nodes.back();
                nodes.pop_back();
                INDEX currentDistance = distances[currentNode];
                if (inverseLocalSpheres[currentDistance].first == -1 && connectedComponentIndices[currentNode] != -1) {
                    inverseLocalSpheres[currentDistance].first = connectedComponentIndices[currentNode];
                }
                if (!visited[currentNode]) {
                    visited[currentNode] = true;
                    inverseLocalSpheres[currentDistance].second.emplace_back(currentNode);
                    for (NodeId i = 0; i < _graph.degree(currentNode); ++i) {
                        NodeId neighbor = _graph.neighbor(currentNode, i);
                        INDEX neighborDistance = distances[neighbor];
                        if (neighborDistance <= currentDistance) {
                            nodes.push_front(neighbor);
                        }
                    }

                }
            }
            // transform the inverse local spheres into a partition
            int counter = 0;
            for (auto & inverseLocalSphere : inverseLocalSpheres) {
                if (inverseLocalSphere.first == -1) {
                    partition[counter].emplace_back(inverseLocalSphere.second);
                    // set the cluster indices
                    for (auto & node : inverseLocalSphere.second) {
                        connectedComponentIndices[node] = (int) partition[counter].size() - 1;
                    }
                }
                else{
                    partition[counter][inverseLocalSphere.first].insert(partition[counter][inverseLocalSphere.first].end(), inverseLocalSphere.second.begin(), inverseLocalSphere.second.end());
                    // set the cluster indices
                    for (auto & node : inverseLocalSphere.second) {
                        connectedComponentIndices[node] = inverseLocalSphere.first;
                    }
                }
                ++counter;


            }
        }
    }





    // create leveling of the graph with respect to the root node
    LayeringPartition layeringPartition;
    // iterate over distance vector
    for (INDEX i = 0; i < distances.size(); ++i) {
        // get the distance
        INDEX distance = distances[i];
        while (layeringPartition.layeringPartition.size() <= distance){
            layeringPartition.layeringPartition.emplace_back();
            layeringPartition.layeringPartition[distance].layeringCluster.emplace_back();
        }
        layeringPartition.layeringPartition[distance].layeringCluster[0].push_back(i);

    }
    return layeringPartition;

}

#endif //LIBGRAPH_TREEEMBEDDING_H
