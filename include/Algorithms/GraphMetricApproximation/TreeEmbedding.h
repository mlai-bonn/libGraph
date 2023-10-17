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

inline LayeringPartition TreeEmbeddingAlgorithm::GetLayeringPartition(const TreeEmbeddingInputParameters& inputParameters){
    std::vector<INDEX> distances;
    GraphStruct::BFSDistances(_graph, inputParameters.root_node_id, distances);

    std::vector<std::vector<NodeId>> spheres;

    std::vector<std::vector<std::vector<NodeId>>> partition;
    std::vector<int> clusterIndices = std::vector<int>(_graph.nodes(), -1);
    for (INDEX i = 0; i < distances.size(); ++i) {
        // get the distance
        INDEX distance = distances[i];
        while (spheres.size() <= distance){
            spheres.emplace_back();
            partition.emplace_back();
        }
        spheres[distance].push_back(i);
    }

    // iterate over the spheres
    for (int level = (int) spheres.size() - 1; level >= 0; --level) {
        auto & sphere = spheres[level];
        // get connected components of the sphere
        std::vector<std::vector<NodeId>> connectedComponents;
        // bfs on the sphere
        std::vector<bool> visited(sphere.size(), false);
        for (auto startNode : sphere) {
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
                    INDEX currentDistance = distances[currentNode];
                    for (NodeId i = 0; i < _graph.degree(currentNode); ++i) {
                        NodeId neighbor = _graph.neighbor(currentNode, i);
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
        for (auto & connectedComponent : connectedComponents) {
            partition[level].emplace_back(connectedComponent);
        }
        // TODO contract the graph
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
                if (inverseLocalSpheres[currentDistance].first == -1 && clusterIndices[currentNode] != -1) {
                    inverseLocalSpheres[currentDistance].first = clusterIndices[currentNode];
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
                        clusterIndices[node] = (int) partition[counter].size() - 1;
                    }
                }
                else{
                    partition[counter][inverseLocalSphere.first].insert(partition[counter][inverseLocalSphere.first].end(), inverseLocalSphere.second.begin(), inverseLocalSphere.second.end());
                    // set the cluster indices
                    for (auto & node : inverseLocalSphere.second) {
                        clusterIndices[node] = inverseLocalSphere.first;
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
