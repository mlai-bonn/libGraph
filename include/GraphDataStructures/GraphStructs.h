//
// Created by florian on 26.10.23.
//

#ifndef GOOGLE_TESTS_GRAPHSTRUCTS_H
#define GOOGLE_TESTS_GRAPHSTRUCTS_H

#import "typedefs.h"
#import <iostream>


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


#endif //GOOGLE_TESTS_GRAPHSTRUCTS_H
