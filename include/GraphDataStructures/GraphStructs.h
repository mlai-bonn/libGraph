//
// Created by florian on 26.10.23.
//

#ifndef GRAPH_STRUCTS_H
#define GRAPH_STRUCTS_H

#include "typedefs.h"
#include <iostream>

/**
 * Defines the type of the graph (general, tree, outerplanar)
 * Some algorithms can be optimized if the type is known in advance
*/
enum class GraphType{
    GENERAL,
    TREE,
    OUTERPLANAR,
};

inline std::ostream& operator<<(std::ostream& os, const GraphType& graphType) {
    switch (graphType) {
        case GraphType::GENERAL:
            os << "GENERAL";
            break;
        case GraphType::TREE:
            os << "TREE";
            break;
        case GraphType::OUTERPLANAR:
            os << "OUTERPLANAR";
            break;
    }
    return os;
}

/**
 * Defines the format of the graph file (for saving and loading)
*/
enum class GraphFormat{
    BGF,
    BGFS,
    BINARY,
    EDGES,
    PEREGRINE_DATA,
    PEREGRINE_SMALL,
    DIMACS,
    AIDS,
};

inline std::ostream& operator<<(std::ostream& os, const GraphFormat& graphFormat) {
    switch (graphFormat) {
        case GraphFormat::BGF:
            os << "BGF";
            break;
        case GraphFormat::BGFS:
            os << "BGFS";
            break;
        case GraphFormat::BINARY:
            os << "BINARY";
            break;
        case GraphFormat::EDGES:
            os << "EDGES";

        default:
            os << "UNKNOWN";
            break;
    }
    return os;
}

/**
 * Struct for saving the graph to a file
 * graphPath: path to the graph file
 * Name: name of the graph
 * Format: format of the graph file
 * Labeled: if true the graph will be saved with labels
 */
struct SaveParams{
    std::string graphPath;
    std::string Name;
    GraphFormat Format = GraphFormat::BGFS;
    bool Labeled = false;
};



struct NodePair{
    NodePair() = default;
    NodePair(NodeId src, NodeId dst, bool sort = true){
        if (sort){
            source = std::min(src, dst);
            destination = std::max(src, dst);
            return;
        }
        else {
            source = src;
            destination = dst;
        }
    };
    explicit NodePair(const std::pair<NodeId, NodeId>& pair, bool sort = true){
        if (sort){
            source = std::min(pair.first, pair.second);
            destination = std::max(pair.first, pair.second);
            return;
        }
        else {
            source = pair.first;
            destination = pair.second;
        }
    };
    bool operator == (const NodePair& nodePair) const{
        return source == nodePair.source && destination == nodePair.destination;
    }
    bool operator < (const NodePair& nodePair) const{
        return source < nodePair.source || source == nodePair.source && destination < nodePair.destination;
    }
    [[nodiscard]] NodeId first() const{return source;};
    [[nodiscard]] NodeId second() const{return destination;};

    void print() const{
        std::cout << " " << source << "-" << destination << " ";
    }

private:
    NodeId source;
    NodeId destination;
};

// The specialized hash function for `unordered_map` keys
struct hashNodePair
{
    std::size_t operator() (const NodePair &nodePair) const
    {
        std::size_t h1 = std::hash<NodeId>()(nodePair.first());
        std::size_t h2 = std::hash<NodeId>()(nodePair.second());
        return h1 ^ h2;
    }
};


#endif //GRAPH_STRUCTS_H
