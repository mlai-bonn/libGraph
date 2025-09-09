//
// Created by florian on 09.09.25.
//

#ifndef GEDEXAMPLE_GEDSTRUCTS_H
#define GEDEXAMPLE_GEDSTRUCTS_H
#include <iosfwd>
#include <unordered_set>

#include "typedefs.h"
#include "GraphDataStructures/GraphBase.h"


enum class EditType {
    INSERT,
    DELETE,
    RELABEL,
  };

inline std::ostream& operator<<(std::ostream& os, const EditType& type) {
    switch (type) {
        case EditType::INSERT:
            return os << "INSERT";
        case EditType::DELETE:
            return os << "DELETE";
        case EditType::RELABEL:
            return os << "RELABEL";
        default:
            return os;
    }
}


enum class OperationObject {
    NODE,
    EDGE,
  };

inline std::ostream& operator<<(std::ostream& os, const OperationObject& operationObject) {
    switch (operationObject) {
        case OperationObject::NODE:
            return os << "NODE";
        case OperationObject::EDGE:
            return os << "EDGE";
        default:
            return os;
    }
}



struct EditOperation {
    OperationObject operationObject;
    EditType type;
    NodeId node = -1;
    EDGE edge = {-1, -1};

    bool operator==(const EditOperation& other) const {
        return operationObject == other.operationObject &&
               type == other.type &&
               node == other.node &&
               edge == other.edge;
    }

    bool operator<(const EditOperation& other) const {
        if (operationObject != other.operationObject) {
            return operationObject < other.operationObject;
        }
        if (type != other.type) {
            return type < other.type;
        }
        if (node != other.node) {
            return node < other.node;
        }
        return edge < other.edge;
    }

};

class EditOperationHash {
public:
    size_t operator()(const EditOperation& operation) const {
        size_t h1 = std::hash<int>()(static_cast<int>(operation.operationObject));
        size_t h2 = std::hash<int>()(static_cast<int>(operation.type));
        size_t h3 = std::hash<NodeId>()(operation.node);
        size_t h4 = std::hash<INDEX>()(operation.edge.first);
        size_t h5 = std::hash<INDEX>()(operation.edge.second);
        return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4); // combine hashes
    }
};



inline std::ostream& operator<<(std::ostream& os, const EditOperation& operation) {
    switch (operation.operationObject) {
        case OperationObject::NODE:
            os << "NODE " << operation.node;
            break;
        case OperationObject::EDGE:
            os << "EDGE " << operation.edge.first << "--" << operation.edge.second;
            break;
    }
    switch (operation.type) {
        case EditType::INSERT:
            os << "INSERT" << std::endl;
        case EditType::DELETE:
            os << "DELETE" << std::endl;
        case EditType::RELABEL:
            os << "RELABEL" << std::endl;
    }
    return os;
}


struct EditPath {
    double total_costs;
    std::unordered_set<EditOperation, EditOperationHash> edit_operations;
    std::vector<EditOperation> sequence_of_operations;
    std::vector<GraphStruct> edit_path_graphs;

    GraphStruct source_graph;
    GraphStruct target_graph;

    // tmp variables
    std::vector<NodeId> source_to_current = std::vector<NodeId>();
    std::vector<NodeId> target_to_current = std::vector<NodeId>();
    std::pair<Nodes, Nodes> node_mapping = {std::vector<NodeId>(), std::vector<NodeId>()};

    std::unordered_set<EditOperation, EditOperationHash> remaining_operations;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edit_operations;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edge_deletions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edge_insertions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_node_deletions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_node_insertions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edge_relabels;
    std::unordered_set<EditOperation, EditOperationHash> remaining_node_relabels;
};


#endif //GEDEXAMPLE_GEDSTRUCTS_H