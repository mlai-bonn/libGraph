//
// Created by florian on 29.08.25.
//

#ifndef GEDEXAMPLE_GEDFUNCTIONS_H
#define GEDEXAMPLE_GEDFUNCTIONS_H
#include <map>

// for approximated GED computation we use the gedlib library see https://github.com/dbblumenthal/gedlib
// the following functions are only to embedd the library results in our code base, i.e., the evaluation of the result

// main function to run the examples


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


struct GEDResult {
    double distance;
    std::pair<GraphStruct&, GraphStruct&> graphs;
    std::pair<Nodes, Nodes> node_mapping;
    double time;
    void get_edit_operations(std::unordered_set<EditOperation, EditOperationHash>& edit_operations) const;
    void get_edit_path(EditPath& edit_path, int seed = 0) const;

    // edit path extension
    void delete_edge(EditPath& edit_path, const EditOperation& operation) const;

    static void insert_edge(EditPath& edit_path, const EditOperation& operation);
    static void relabel_edge(EditPath& edit_path, const EditOperation& operation);
    void delete_node(EditPath& edit_path, const EditOperation& operation) const;

    static void insert_node(EditPath& edit_path, const EditOperation& operation);

    static void relabel_node(EditPath& edit_path, const EditOperation& operation);
    
};

inline void GEDResult::get_edit_operations(std::unordered_set<EditOperation, EditOperationHash>& edit_operations) const {
    edit_operations.clear();
    // get list of edit operations from the node mapping (get node operations) and derive edge operations
    const Nodes& mapping_first = node_mapping.first;
    const Nodes& mapping_second = node_mapping.second;


    std::set<NodeId> deleted_nodes;
    std::set<NodeId> inserted_nodes;
    std::set<NodeId> relabeled_nodes;
    std::set<EDGE> deleted_edges;
    std::set<EDGE> inserted_edges;
    std::set<EDGE> relabeled_edges;
    // node operations

    // Find deleted nodes
    NodeId source_id = 0;
    for (auto target_id : mapping_first) {
        if (target_id > mapping_second.size()) {
            deleted_nodes.insert(source_id);
        }
        ++source_id;
    }

    // Find inserted nodes
    NodeId target_id = 0;
    for (auto source_id : mapping_second) {
        if (source_id > mapping_first.size()) {
            inserted_nodes.insert(target_id);
        }
        ++target_id;
    }

    // Find relabeled nodes
    if (graphs.first.labelType != LABEL_TYPE::UNLABELED && graphs.second.labelType != LABEL_TYPE::UNLABELED) {
        source_id = 0;
        for (auto x : mapping_first) {
            if (x < mapping_second.size()) {
                // check if labels are different
                if (graphs.first.label(source_id) != graphs.second.label(x)) {
                    relabeled_nodes.insert(source_id);
                }
            }
            ++source_id;
        }
    }

    // Edge operations
    // Find deleted edges
    // iterate over edges in first graph
    for (NodeId source_i = 0; source_i < graphs.first.nodes(); ++source_i) {
        for (NodeId source_j : graphs.first.get_neighbors(source_i)) {
            NodeId target_i = mapping_first[source_i];
            NodeId target_j = mapping_first[source_j];
            if (source_i < source_j) { // to avoid double counting
                // check if edge exists in second graph
                if (target_i < mapping_second.size() && target_j < mapping_second.size()) {
                    if (!graphs.second.edge(target_i, target_j)) {
                        deleted_edges.insert({source_i, source_j});
                    }
                } else {
                    deleted_edges.insert({source_i, source_j});
                }
            }
        }
    }
    // Find inserted edges
    for (auto target_i = 0; target_i < graphs.second.nodes(); ++target_i) {
        for (auto target_j : graphs.second.get_neighbors(target_i)) {
            if (target_i < target_j) {
                NodeId source_i = mapping_second[target_i];
                NodeId source_j = mapping_second[target_j];
                if (source_i < mapping_first.size() && source_j < mapping_first.size()) {
                    if (!graphs.first.edge(source_i, source_j)) {
                        inserted_edges.insert({target_i, target_j});
                    }
                } else {
                    inserted_edges.insert({target_i, target_j});
                }
            }
        }
    }
    // Find relabeled edges TODO
    // Create the edit operations
    for (auto x : deleted_nodes) {
        EditOperation operation = {
            .operationObject = OperationObject::NODE,
            .type = EditType::DELETE,
            .node = x,
        };
        edit_operations.insert(operation);
    }
    for (auto x : inserted_nodes) {
        EditOperation operation = {
            .operationObject = OperationObject::NODE,
            .type = EditType::INSERT,
            .node = x,
        };
        edit_operations.insert(operation);
    }
    for (auto x : relabeled_nodes) {
        EditOperation operation = {
            .operationObject = OperationObject::NODE,
            .type = EditType::RELABEL,
            .node = x,
        };
        edit_operations.insert(operation);
    }
    for (auto x : deleted_edges) {
        EditOperation operation = {
            .operationObject = OperationObject::EDGE,
            .type = EditType::DELETE,
            .edge = x,
        };
        edit_operations.insert(operation);
    }
    for (auto x : inserted_edges) {
        EditOperation operation = {
            .operationObject = OperationObject::EDGE,
            .type = EditType::INSERT,
            .edge = x,
        };
        edit_operations.insert(operation);
    }
    for (auto x : relabeled_edges) {
        EditOperation operation = {
            .operationObject = OperationObject::EDGE,
            .type = EditType::RELABEL,
            .edge = x,
        };
        edit_operations.insert(operation);
    }
}

inline void GEDResult::get_edit_path(EditPath& edit_path, int seed) const {
    edit_path.total_costs = distance;
    edit_path.edit_operations.clear();
    get_edit_operations(edit_path.edit_operations);


    // add first graph to edit path graphs
    edit_path.source_graph = graphs.first;
    edit_path.target_graph = graphs.second;

    edit_path.edit_path_graphs.clear();
    edit_path.edit_path_graphs.push_back(graphs.first);
    edit_path.source_to_current = std::vector<NodeId>(graphs.first.nodes(), 0);
    std::iota(std::begin(edit_path.source_to_current), std::end(edit_path.source_to_current), 0);
    edit_path.target_to_current = node_mapping.second;
    edit_path.node_mapping = node_mapping;

    edit_path.remaining_operations = edit_path.edit_operations;
    

    for (auto x : edit_path.remaining_operations) {
        if (x.operationObject == OperationObject::NODE) {
            if (x.type == EditType::DELETE) {
                edit_path.remaining_node_deletions.insert(x);
            } else if (x.type == EditType::INSERT) {
                edit_path.remaining_node_insertions.insert(x);
            } else if (x.type == EditType::RELABEL) {
                edit_path.remaining_node_relabels.insert(x);
            }
        } else if (x.operationObject == OperationObject::EDGE) {
            if (x.type == EditType::DELETE) {
                edit_path.remaining_edge_deletions.insert(x);
            } else if (x.type == EditType::INSERT) {
                edit_path.remaining_edge_insertions.insert(x);
            } else if (x.type == EditType::RELABEL) {
                edit_path.remaining_edge_relabels.insert(x);
            }
        }
    }

    // Get a meaningful sequence of edit_path operations
    while (!edit_path.remaining_operations.empty()) {
        // first edge deletions
        while (!edit_path.remaining_edge_deletions.empty()) {
            auto it = edit_path.remaining_edge_deletions.begin();
            delete_edge(edit_path, *it);
        }
        // Then edge insertions
        while (!edit_path.remaining_edge_insertions.empty()) {
            auto it = edit_path.remaining_edge_insertions.begin();
            insert_edge(edit_path, *it);
        }

        // Then remaining node deletions
        while (!edit_path.remaining_node_deletions.empty()) {
            auto it = edit_path.remaining_node_deletions.begin();
            delete_node(edit_path, *it);
        }

        // Then remaining node insertions
        while (!edit_path.remaining_node_insertions.empty()) {
            auto it = edit_path.remaining_node_insertions.begin();
            insert_node(edit_path, *it);
        }
        // TODO relabellings
        while (!edit_path.remaining_node_relabels.empty()) {
            auto it = edit_path.remaining_node_relabels.begin();
            relabel_node(edit_path, *it);
        }
        while (!edit_path.remaining_edge_relabels.empty()) {
            auto it = edit_path.remaining_edge_relabels.begin();
            relabel_edge(edit_path, *it);
        }

    }


}

inline void GEDResult::delete_edge(EditPath &edit_path, const EditOperation &operation) const {
    const NodeId source_i = operation.edge.first;
    const NodeId source_j = operation.edge.second;
    const NodeId current_i = edit_path.source_to_current[source_i];
    const NodeId current_j = edit_path.source_to_current[source_j];
    // if an edge deletion creates an isolated node, delete the node as well
    // new graph after deleting the edge
    GraphStruct new_graph = edit_path.edit_path_graphs.back();
    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);
    new_graph.remove_edge(current_i, current_j);
    edit_path.edit_path_graphs.emplace_back(new_graph);

    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_operations.erase(operation);
    edit_path.remaining_edge_deletions.erase(operation);
    // check if node1 is isolated
    if (new_graph.degree(current_i) == 0) {
        for (const auto x : edit_path.remaining_node_deletions) {
            if (x.node == source_i) {
                delete_node(edit_path, x);
                break;
            }
        }
    }
    // check if node2 is isolated
    if (new_graph.degree(current_j) == 0) {
        for (const auto x : edit_path.remaining_node_deletions) {
            if (x.node == source_j) {
                delete_node(edit_path, x);
                break;
            }
        }
    }
}

inline void GEDResult::delete_node(EditPath &edit_path, const EditOperation &operation) const {
    const NodeId source_node = operation.node;
    const NodeId current_node = edit_path.source_to_current[source_node];
    // TODO check whether node has degree 0 o.w. delete all edges
    if (edit_path.edit_path_graphs.back().degree(current_node) != 0) {
        // raise an error
        std::cerr << "Error: Trying to delete a node that is not isolated. Node: " << source_node << " Current Target Node: " << current_node << std::endl;

    }
    // new graph after deleting the node
    GraphStruct new_graph;
    // build up new graph from last graph without the node
    new_graph = GraphStruct(edit_path.edit_path_graphs.back().nodes() - 1, Labels());
    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);

    // Updating the maps
    for (NodeId i = source_node + 1; i < graphs.first.nodes(); ++i) {
        edit_path.source_to_current[i] -= 1;
    }
    for (NodeId j = 0; j < graphs.second.nodes(); ++j) {
        if (edit_path.target_to_current[j] > current_node) {
            edit_path.target_to_current[j] -=1;
        }
    }
    edit_path.source_to_current[source_node] = -1;

    // add the edges between the correct nodes
    for (INDEX i = 0; i < edit_path.edit_path_graphs.back().nodes(); ++i) {
        for (const auto j : edit_path.edit_path_graphs.back().get_neighbors(i)) {
            if (i < j) {
                if (i>= current_node) {
                    new_graph.add_edge(i - 1, j-1);
                }
                else {
                    if (j >= current_node) {
                        new_graph.add_edge(i, j - 1);
                    }
                    else {
                        new_graph.add_edge(i, j);
                    }
                }
            }
        }
    }
    // Add the labels
    if (edit_path.edit_path_graphs.back().labelType != LABEL_TYPE::UNLABELED) {
        Labels labels = std::vector<Label>(edit_path.edit_path_graphs.back().labels().size() - 1, 0);
        for (NodeId i = 0, j = 0; i < edit_path.edit_path_graphs.back().labels().size(); ++i) {
            if (i != current_node) {
                labels[j] = edit_path.edit_path_graphs.back().label(i);
                ++j;
            }
        }
        new_graph.set_labels(&labels);
    }

    edit_path.edit_path_graphs.emplace_back(new_graph);
    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_node_deletions.erase(operation);
    edit_path.remaining_operations.erase(operation);
}


inline void GEDResult::insert_edge(EditPath &edit_path, const EditOperation &operation) {
    // Check if both nodes of the edge are already inserted
    const NodeId target_i = operation.edge.first;
    const NodeId target_j = operation.edge.second;

    NodeId current_i = edit_path.target_to_current[target_i];
    NodeId current_j = edit_path.target_to_current[target_j];

    const bool first_node_inserted = current_i <= edit_path.edit_path_graphs.back().nodes();
    const bool second_node_inserted = current_j <= edit_path.edit_path_graphs.back().nodes();

    if (!first_node_inserted) {
        const EditOperation first_node_insertion = {
            .operationObject = OperationObject::NODE,
            .type = EditType::INSERT,
            .node = operation.edge.first,
        };
        insert_node(edit_path, first_node_insertion);
    }
    if (!second_node_inserted) {
        const EditOperation second_node_insertion = {
            .operationObject = OperationObject::NODE,
            .type = EditType::INSERT,
            .node = operation.edge.second,
        };
        insert_node(edit_path, second_node_insertion);
    }
    // new graph after inserting the edge
    current_i = edit_path.target_to_current[target_i];
    current_j = edit_path.target_to_current[target_j];
    GraphStruct new_graph = edit_path.edit_path_graphs.back();
    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);

    new_graph.add_edge(current_i, current_j);
    edit_path.edit_path_graphs.emplace_back(new_graph);


    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_operations.erase(operation);
    edit_path.remaining_edge_insertions.erase(operation);
}

inline void GEDResult::insert_node(EditPath &edit_path, const EditOperation &operation) {
    // new graph after inserting the node
    auto new_graph = GraphStruct(edit_path.edit_path_graphs.back().nodes() +1, Labels());

    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);


    //Update the maps
    edit_path.target_to_current[operation.node] = edit_path.edit_path_graphs.back().nodes();

    for (INDEX i = 0; i < edit_path.edit_path_graphs.back().nodes(); ++i) {
        for (auto j : edit_path.edit_path_graphs.back().get_neighbors(i)) {
            if (i < j) {
                new_graph.add_edge(i, j);
            }
        }
    }

    if (edit_path.edit_path_graphs.back().labelType != LABEL_TYPE::UNLABELED) {
        Labels labels = edit_path.edit_path_graphs.back().labels();
        labels.push_back(edit_path.target_graph.label(operation.node));
        new_graph.set_labels(&labels);
    }

    edit_path.edit_path_graphs.emplace_back(new_graph);



    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_node_insertions.erase(operation);
    edit_path.remaining_operations.erase(operation);
}

inline void GEDResult::relabel_edge(EditPath &edit_path, const EditOperation &operation) {
    // TODO
    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_edge_relabels.erase(operation);
    edit_path.remaining_operations.erase(operation);
}

inline void GEDResult::relabel_node(EditPath &edit_path, const EditOperation &operation) {
    const NodeId source_i = operation.node;
    const NodeId current_i = edit_path.source_to_current[source_i];
    const NodeId target_i = edit_path.node_mapping.first[source_i];
    Label target_label = edit_path.target_graph.label(target_i);
    GraphStruct new_graph = edit_path.edit_path_graphs.back();
    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);
    // relabel the node
    Labels labels = new_graph.labels();
    labels[current_i] = target_label;
    new_graph.set_labels(&labels);
    edit_path.edit_path_graphs.emplace_back(new_graph);

    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_operations.erase(operation);
    edit_path.remaining_node_relabels.erase(operation);
}



#endif //GEDEXAMPLE_GEDFUNCTIONS_H