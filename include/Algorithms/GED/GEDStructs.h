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
    DeleteEdges,
    InsertEdges,
    DeleteNodes,
    InsertNodes,
    RelabelEdges,
    RelabelNodes,
    RandomDeleteEdges,
    RandomInsertEdges,
    RandomDeleteNodes,
    RandomInsertNodes,
    RandomRelabelEdges,
    RandomRelabelNodes,
    DeleteIsolatedNodes,
  };

inline std::string EditPathStrategyToString(const EditPathStrategy& strategy) {
    switch (strategy) {
        case EditPathStrategy::Random:
            return "Random";
        case EditPathStrategy::DeleteEdges:
            return "DeleteEdges";
        case EditPathStrategy::InsertEdges:
            return "InsertEdges";
        case EditPathStrategy::DeleteNodes:
            return "DeleteNodes";
        case EditPathStrategy::InsertNodes:
            return "InsertNodes";
        case EditPathStrategy::RelabelEdges:
            return "RelabelEdges";
        case EditPathStrategy::RelabelNodes:
            return "RelabelNodes";
        case EditPathStrategy::RandomDeleteEdges:
            return "RandomDeleteEdges";
        case EditPathStrategy::RandomInsertEdges:
            return "RandomInsertEdges";
        case EditPathStrategy::RandomDeleteNodes:
            return "RandomDeleteNodes";
        case EditPathStrategy::RandomInsertNodes:
            return "RandomInsertNodes";
        case EditPathStrategy::RandomRelabelEdges:
            return "RandomRelabelEdges";
        case EditPathStrategy::RandomRelabelNodes:
            return "RandomRelabelNodes";
        case EditPathStrategy::DeleteIsolatedNodes:
            return "DeleteIsolatedNodes";
        default:
            return "Random";
    }
}

inline std::string EditPathStrategyToStringShort(const EditPathStrategy& strategy) {
    switch (strategy) {
        case EditPathStrategy::Random:
            return "Rnd";
        case EditPathStrategy::DeleteEdges:
            return "d-E";
        case EditPathStrategy::InsertEdges:
            return "i-E";
        case EditPathStrategy::DeleteNodes:
            return "d-N";
        case EditPathStrategy::InsertNodes:
            return "i-N";
        case EditPathStrategy::RelabelEdges:
            return "r-E";
        case EditPathStrategy::RelabelNodes:
            return "r-N";
        case EditPathStrategy::RandomDeleteEdges:
            return "Rnd-d-E";
        case EditPathStrategy::RandomInsertEdges:
            return "Rnd-i-E";
        case EditPathStrategy::RandomDeleteNodes:
            return "Rnd-d-N";
        case EditPathStrategy::RandomInsertNodes:
            return "Rnd-i-N";
        case EditPathStrategy::RandomRelabelEdges:
            return "Rnd-r-E";
        case EditPathStrategy::RandomRelabelNodes:
            return "Rnd-r-N";
        case EditPathStrategy::DeleteIsolatedNodes:
            return "d-IsoN";
        default:
            return "Rnd";
    }
}

inline std::string EditPathStrategiesToStringShort(const std::vector<EditPathStrategy>& strategies) {
    // For naming output folders
    std::string result;
    // join strategies with _
    for (size_t i = 0; i < strategies.size(); ++i) {
        result += EditPathStrategyToStringShort(strategies[i]);
        if (i < strategies.size() - 1) {
            result += "_";
        }
    }
    return result;
}

inline std::vector<std::string> EditPathStrategiesToString(const std::vector<EditPathStrategy>& strategies) {
    std::vector<std::string> strategy_strings;
    for (const auto& strategy : strategies) {
        strategy_strings.push_back(EditPathStrategyToString(strategy));
    }
    return strategy_strings;
}

inline EditPathStrategy StringToEditPathStrategy(const std::string& strategy_str) {
    if (strategy_str == "Random") {
        return EditPathStrategy::Random;
    } else if (strategy_str == "DeleteEdges") {
        return EditPathStrategy::DeleteEdges;
    } else if (strategy_str == "InsertEdges") {
        return EditPathStrategy::InsertEdges;
    } else if (strategy_str == "DeleteNodes") {
        return EditPathStrategy::DeleteNodes;
    } else if (strategy_str == "InsertNodes") {
        return EditPathStrategy::InsertNodes;
    } else if (strategy_str == "DeleteIsolatedNodes") {
        return EditPathStrategy::DeleteIsolatedNodes;
    } else if (strategy_str == "RelabelEdges") {
        return EditPathStrategy::RelabelEdges;
    } else if (strategy_str == "RelabelNodes") {
        return EditPathStrategy::RelabelNodes;
    } else if (strategy_str == "RandomDeleteEdges") {
        return EditPathStrategy::RandomDeleteEdges;
    } else if (strategy_str == "RandomInsertEdges") {
        return EditPathStrategy::RandomInsertEdges;
    } else if (strategy_str == "RandomDeleteNodes") {
        return EditPathStrategy::RandomDeleteNodes;
    } else if (strategy_str == "RandomInsertNodes") {
        return EditPathStrategy::RandomInsertNodes;
    } else if (strategy_str == "RandomRelabelEdges") {
        return EditPathStrategy::RandomRelabelEdges;
    } else if (strategy_str == "RandomRelabelNodes") {
        return EditPathStrategy::RandomRelabelNodes;
    } else {
        // print error and list valid strategies
        std::cerr << "Error: Invalid edit path strategy: " << strategy_str << std::endl;
        std::cerr << "Valid strategies are: Random, DeleteEdges, InsertEdges, DeleteNodes, InsertNodes, RelabelEdges, RelabelNodes, RandomDeleteEdges, RandomInsertEdges, RandomDeleteNodes, RandomInsertNodes, RandomRelabelEdges, RandomRelabelNodes, DeleteIsolatedNodes" << std::endl;
        return EditPathStrategy::Random;
    }
}

inline std::vector<EditPathStrategy> StringsToEditPathStrategies(const std::vector<std::string>& strategy_strs) {
    std::vector<EditPathStrategy> strategies;
    for (const auto& strategy_str : strategy_strs) {
        strategies.push_back(StringToEditPathStrategy(strategy_str));
    }
    return strategies;
}

inline std::vector<EditPathStrategy> StringToEditPathStrategies(const std::string& strategies_str) {
    // split by space, comma or semicolon
    std::vector<std::string> strategy_strs;
    size_t start = 0;
    size_t end = strategies_str.find_first_of(" ,;");
    while (end != std::string::npos) {
        strategy_strs.push_back(strategies_str.substr(start, end - start));
        start = end + 1;
        end = strategies_str.find_first_of(" ,;", start);
    }
    strategy_strs.push_back(strategies_str.substr(start));
    return StringsToEditPathStrategies(strategy_strs);
}



inline bool GetValidStrategy(std::vector<EditPathStrategy>& strategies) {
    // must contain only Random (Random with DeleteIsolatedNodes) or some combination of the six operations (each at most once) if an operation is missing add Random at the end
    bool has_random = false;
    bool has_delete_isolated_nodes = false;

    bool has_delete_edges = false;
    bool has_insert_edges = false;
    bool has_delete_nodes = false;
    bool has_insert_nodes = false;
    bool has_relabel_edges = false;
    bool has_relabel_nodes = false;
    bool has_random_delete_edges = false;
    bool has_random_insert_edges = false;
    bool has_random_delete_nodes = false;
    bool has_random_insert_nodes = false;
    bool has_random_relabel_edges = false;
    bool has_random_relabel_nodes = false;

    int counter = 0;
    while (counter < strategies.size()) {
        auto strategy = strategies[counter];
        ++counter;
        switch (strategy) {
            case EditPathStrategy::Random:
                if (has_random) {
                    // remove duplicate Random strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                has_random = true;
                if (strategies.size() == 2) {
                    if (strategies[0] != EditPathStrategy::DeleteIsolatedNodes && strategies[1] != EditPathStrategy::DeleteIsolatedNodes && strategies[0] != EditPathStrategy::Random && strategies[1] != EditPathStrategy::Random) {
                        std::cerr << "Error: Random strategy can only be combined with DeleteIsolatedNodes." << std::endl;
                        return false;
                    }
                }
                if (strategies.size() > 2) {
                    std::cerr << "Error: Random strategy can only be combined with DeleteIsolatedNodes." << std::endl;
                    return false;
                }
                break;
            case EditPathStrategy::DeleteEdges:
                if (has_delete_edges) {
                    // remove duplicate DeleteEdges strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_random_delete_edges) {
                    std::cerr << "Error: DeleteEdges strategy cannot be combined with RandomDeleteEdges." << std::endl;
                    return false;
                }
                has_delete_edges = true;
                break;
            case EditPathStrategy::InsertEdges:
                if (has_insert_edges) {
                    // remove duplicate InsertEdges strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_random_insert_edges) {
                    std::cerr << "Error: InsertEdges strategy cannot be combined with RandomInsertEdges." << std::endl;
                    return false;
                }
                has_insert_edges = true;
                break;
            case EditPathStrategy::DeleteNodes:
                if (has_delete_nodes) {
                    // remove duplicate DeleteNodes strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_random_delete_nodes) {
                    std::cerr << "Error: DeleteNodes strategy cannot be combined with RandomDeleteNodes." << std::endl;
                    return false;
                }
                has_delete_nodes = true;
                break;
            case EditPathStrategy::InsertNodes:
                if (has_insert_nodes) {
                    // remove duplicate InsertNodes strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_random_insert_nodes) {
                    std::cerr << "Error: InsertNodes strategy cannot be combined with RandomInsertNodes." << std::endl;
                    return false;
                }
                has_insert_nodes = true;
                break;
            case EditPathStrategy::RelabelEdges:
                if (has_relabel_edges) {
                    // remove duplicate RelabelEdges strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_random_relabel_edges) {
                    std::cerr << "Error: RelabelEdges strategy cannot be combined with RandomRelabelEdges." << std::endl;
                    return false;
                }
                has_relabel_edges = true;
                break;
            case EditPathStrategy::RelabelNodes:
                if (has_relabel_nodes) {
                    // remove duplicate RelabelNodes strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_random_relabel_nodes) {
                    std::cerr << "Error: RelabelNodes strategy cannot be combined with RandomRelabelNodes." << std::endl;
                    return false;
                }
                has_relabel_nodes = true;
                break;
            case EditPathStrategy::RandomDeleteEdges:
                if (has_random_delete_edges) {
                    // remove duplicate RandomDeleteEdges strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_delete_edges) {
                    std::cerr << "Error: RandomDeleteEdges strategy cannot be combined with DeleteEdges." << std::endl;
                    return false;
                }
                has_random_delete_edges = true;
                break;
            case EditPathStrategy::RandomInsertEdges:
                if (has_random_insert_edges) {
                    // remove duplicate RandomInsertEdges strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_insert_edges) {
                    std::cerr << "Error: RandomInsertEdges strategy cannot be combined with InsertEdges." << std::endl;
                    return false;
                }
                has_random_insert_edges = true;
                break;
            case EditPathStrategy::RandomDeleteNodes:
                if (has_random_delete_nodes) {
                    // remove duplicate RandomDeleteNodes strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_delete_nodes) {
                    std::cerr << "Error: RandomDeleteNodes strategy cannot be combined with DeleteNodes." << std::endl;
                    return false;
                }
                has_random_delete_nodes = true;
                break;
            case EditPathStrategy::RandomInsertNodes:
                if (has_random_insert_nodes) {
                   // remove duplicate RandomInsertNodes strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_insert_nodes) {
                    std::cerr << "Error: RandomInsertNodes strategy cannot be combined with InsertNodes." << std::endl;
                    return false;
                }
                has_random_insert_nodes = true;
                break;
            case EditPathStrategy::RandomRelabelEdges:
                if (has_random_relabel_edges) {
                    // remove duplicate RandomRelabelEdges strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_relabel_edges) {
                    std::cerr << "Error: RandomRelabelEdges strategy cannot be combined with RelabelEdges." << std::endl;
                    return false;
                }
                has_random_relabel_edges = true;
                break;
            case EditPathStrategy::RandomRelabelNodes:
                if (has_random_relabel_nodes) {
                    // remove duplicate RandomRelabelNodes strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                if (has_relabel_nodes) {
                    std::cerr << "Error: RandomRelabelNodes strategy cannot be combined with RelabelNodes." << std::endl;
                    return false;
                }
                has_random_relabel_nodes = true;
                break;
            case EditPathStrategy::DeleteIsolatedNodes:
                if (has_delete_isolated_nodes) {
                    // remove duplicate DeleteIsolatedNodes strategies
                    strategies.erase(strategies.begin() + counter - 1);
                    --counter;
                    continue;
                }
                has_delete_isolated_nodes = true;
                break;
        }
    }
    // If empty add Random
    if (strategies.empty()) {
        strategies.push_back(EditPathStrategy::Random);
        return true;
    }
    if (has_random) {
            return  true;
    }
     // If some strategies are not present the default is their random version
    if (!has_delete_edges) {
        strategies.push_back(EditPathStrategy::RandomDeleteEdges);
    }
    if (!has_insert_edges) {
        strategies.push_back(EditPathStrategy::RandomInsertEdges);
    }
    if (!has_delete_nodes) {
        strategies.push_back(EditPathStrategy::RandomDeleteNodes);
    }
    if (!has_insert_nodes) {
        strategies.push_back(EditPathStrategy::RandomInsertNodes);
    }
    if (!has_relabel_edges) {
        strategies.push_back(EditPathStrategy::RandomRelabelEdges);
    }
    if (!has_relabel_nodes) {
        strategies.push_back(EditPathStrategy::RandomRelabelNodes);
    }
    return true;
}



inline std::ostream& operator<<(std::ostream& os, const EditPathStrategy& strategy) {
    switch (strategy) {
        case EditPathStrategy::Random:
            return os << "Random";
        case EditPathStrategy::DeleteEdges:
            return os << "DeleteEdges";
        case EditPathStrategy::InsertEdges:
            return os << "InsertEdges";
        case EditPathStrategy::DeleteNodes:
            return os << "DeleteNodes";
        case EditPathStrategy::InsertNodes:
            return os << "InsertNodes";
        case EditPathStrategy::RelabelEdges:
            return os << "RelabelEdges";
        case EditPathStrategy::RelabelNodes:
            return os << "RelabelNodes";
        case EditPathStrategy::RandomDeleteEdges:
            return os << "RandomDeleteEdges";
        case EditPathStrategy::RandomInsertEdges:
            return os << "RandomInsertEdges";
        case EditPathStrategy::RandomDeleteNodes:
            return os << "RandomDeleteNodes";
        case EditPathStrategy::RandomInsertNodes:
            return os << "RandomInsertNodes";
        case EditPathStrategy::RandomRelabelEdges:
            return os << "RandomRelabelEdges";
        case EditPathStrategy::RandomRelabelNodes:
            return os << "RandomRelabelNodes";
        case EditPathStrategy::DeleteIsolatedNodes:
            return os << "DeleteIsolatedNodes";
        default:
            return os << "Random";
    }
}

inline std::ostream& operator<<(std::ostream& os, const std::vector<EditPathStrategy>& strategies) {
    for (size_t i = 0; i < strategies.size(); ++i) {
        os << strategies[i];
        if (i < strategies.size() - 1) {
            os << ", ";
        }
    }
    return os;
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