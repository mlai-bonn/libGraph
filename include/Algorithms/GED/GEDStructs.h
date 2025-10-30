//
// Created by florian on 09.09.25.
//

#ifndef GED_STRUCTS_H
#define GED_STRUCTS_H
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

enum class EditPathStrategy {
    Random,
    DeleteFirstEdgesNodes,
    InsertFirstEdgesNodes,
    DeleteFirstNodesEdges,
    InsertFirstNodesEdges,
  };

inline std::string EditPathStrategyToString(const EditPathStrategy& strategy) {
    switch (strategy) {
        case EditPathStrategy::Random:
            return "Random";
        case EditPathStrategy::DeleteFirstEdgesNodes:
            return "DeleteFirstEdgesNodes";
        case EditPathStrategy::InsertFirstEdgesNodes:
            return "InsertFirstEdgesNodes";
        case EditPathStrategy::DeleteFirstNodesEdges:
            return "DeleteFirstNodesEdges";
        case EditPathStrategy::InsertFirstNodesEdges:
            return "InsertFirstNodesEdges";
        default:
            return "Unknown";
    }
}

inline EditPathStrategy StringToEditPathStrategy(const std::string& strategy_str) {
    if (strategy_str == "Random") {
        return EditPathStrategy::Random;
    } else if (strategy_str == "DeleteFirstEdgesNodes") {
        return EditPathStrategy::DeleteFirstEdgesNodes;
    } else if (strategy_str == "InsertFirstEdgesNodes") {
        return EditPathStrategy::InsertFirstEdgesNodes;
    } else if (strategy_str == "DeleteFirstNodesEdges") {
        return EditPathStrategy::DeleteFirstNodesEdges;
    } else if (strategy_str == "InsertFirstNodesEdges") {
        return EditPathStrategy::InsertFirstNodesEdges;
    } else {
        throw std::invalid_argument("Unknown EditPathStrategy string: " + strategy_str);
    }
}

inline std::ostream& operator<<(std::ostream& os, const EditPathStrategy& strategy) {
    switch (strategy) {
        case EditPathStrategy::Random:
            return os << "Random";
        case EditPathStrategy::DeleteFirstEdgesNodes:
            return os << "DeleteFirstEdgesNodes";
        case EditPathStrategy::InsertFirstEdgesNodes:
            return os << "InsertFirstEdgesNodes";
        case EditPathStrategy::DeleteFirstNodesEdges:
            return os << "DeleteFirstNodesEdges";
        case EditPathStrategy::InsertFirstNodesEdges:
            return os << "InsertFirstNodesEdges";
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

    void WriteToBinary(std::ofstream& file) const {
        file.write(reinterpret_cast<const char *>(&operationObject), sizeof(operationObject));
        file.write(reinterpret_cast<const char *>(&type), sizeof(type));
        file.write(reinterpret_cast<const char *>(&node), sizeof(node));
        file.write(reinterpret_cast<const char *>(&edge.first), sizeof(edge.first));
        file.write(reinterpret_cast<const char *>(&edge.second), sizeof(edge.second));
    }

    void ReadFromBinary(std::ifstream& file) {
        file.read(reinterpret_cast<char *>(&operationObject), sizeof(operationObject));
        file.read(reinterpret_cast<char *>(&type), sizeof(type));
        file.read(reinterpret_cast<char *>(&node), sizeof(node));
        file.read(reinterpret_cast<char *>(&edge.first), sizeof(edge.first));
        file.read(reinterpret_cast<char *>(&edge.second), sizeof(edge.second));
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
        default:
            // raise error
            os << "UNKNOWN OBJECT";
            break;
    }
    switch (operation.type) {
        case EditType::INSERT:
            os << " INSERT" << std::endl;
            break;
        case EditType::DELETE:
            os << " DELETE" << std::endl;
            break;
        case EditType::RELABEL:
            os << " RELABEL" << std::endl;
            break;
    }
    return os;
}

template<typename T>
struct EditPath {
    double total_costs;
    std::unordered_set<EditOperation, EditOperationHash> edit_operations;
    std::vector<EditOperation> sequence_of_operations;
    std::vector<EditOperation> sequence_of_operations_current;
    std::vector<T> edit_path_graphs;

    T source_graph;
    T target_graph;

    // tmp variables
    std::vector<NodeId> source_to_current = std::vector<NodeId>();
    std::vector<NodeId> target_to_current = std::vector<NodeId>();
    std::pair<Nodes, Nodes> node_mapping = {std::vector<NodeId>(), std::vector<NodeId>()};

    std::list<EditOperation> remaining_operations;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edit_operations;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edge_deletions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edge_insertions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_node_deletions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_node_insertions;
    std::unordered_set<EditOperation, EditOperationHash> remaining_edge_relabels;
    std::unordered_set<EditOperation, EditOperationHash> remaining_node_relabels;

    void Update(T& graph, const EditOperation& operation, const EditOperation& operation_current);
};

// print Edit Path using the sequence of operations
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const EditPath<T>& edit_path) {
    os << "Edit Costs: " << edit_path.total_costs << std::endl;
    os << "Sequence of Operations: " << std::endl;
    for (const auto& operation : edit_path.sequence_of_operations) {
        os << operation;
    }
    return os;
}

template<typename T>
void EditPath<T>::Update(T& graph, const EditOperation& operation, const EditOperation& operation_current) {
    this->edit_path_graphs.emplace_back(graph);
    this->sequence_of_operations.push_back(operation);
    this->sequence_of_operations_current.push_back(operation_current);
    this->remaining_operations.remove(operation);
    if (operation.operationObject == OperationObject::NODE) {
        if (operation.type == EditType::DELETE) {
            this->remaining_node_deletions.erase(operation);
        } else if (operation.type == EditType::INSERT) {
            this->remaining_node_insertions.erase(operation);
        } else if (operation.type == EditType::RELABEL) {
            this->remaining_node_relabels.erase(operation);
        }
    } else if (operation.operationObject == OperationObject::EDGE) {
        if (operation.type == EditType::DELETE) {
            this->remaining_edge_deletions.erase(operation);
        } else if (operation.type == EditType::INSERT) {
            this->remaining_edge_insertions.erase(operation);
        } else if (operation.type == EditType::RELABEL) {
            this->remaining_edge_relabels.erase(operation);
        }
    }
}

#endif //GED_STRUCTS_H