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
    for (INDEX i = 0; i < distances.size(); ++i) {
        // get the distance
        INDEX distance = distances[i];
        while (spheres.size() <= distance){
            spheres.emplace_back();
            partition.emplace_back();
        }
        spheres[distance].push_back(i);

    }

    std::vector<bool> visited(_graph.nodes(), false);
    // start a backward search from the most distant nodes
    std::vector<NodeId> startNodes = spheres[spheres.size() - 1];
    INDEX clusterIndex = 0;
    while (!startNodes.empty()) {
        for (auto & x : partition){
            x.emplace_back();
        }
        auto startNode = startNodes.back();
        startNodes.pop_back();
        std::deque<NodeId> nodes = std::deque<NodeId>();
        nodes.push_back(startNode);
        while (!nodes.empty()) {
            NodeId currentNode = nodes.back();
            nodes.pop_back();
            INDEX currentDistance = distances[currentNode];
            if (!visited[currentNode]) {
                visited[currentNode] = true;
                partition[currentDistance][clusterIndex].emplace_back(currentNode);

                for (NodeId i = 0; i < _graph.degree(currentNode); ++i) {
                    NodeId neighbor = _graph.neighbor(currentNode, i);
                    if (!visited[neighbor]) {
                        INDEX neighborDistance = distances[neighbor];
                        if (neighborDistance <= currentDistance) {
                            nodes.push_front(neighbor);
                        }
                    }
                }
            }
        }
        ++clusterIndex;
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
