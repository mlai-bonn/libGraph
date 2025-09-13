//
// Created by florian on 09.09.25.
//

#ifndef GEDEXAMPLE_GEDEVALUATION_H
#define GEDEXAMPLE_GEDEVALUATION_H
#include <utility>
#include "GEDStructs.h"
#include "typedefs.h"
#include "GraphDataStructures/GraphBase.h"

template <typename T>
struct GEDEvaluation {
    double distance = std::numeric_limits<double>::max();
    double lower_bound = std::numeric_limits<double>::min();
    double upper_bound = std::numeric_limits<double>::max();
    std::pair<T, T> graphs = {T(), T()};
    std::pair<INDEX, INDEX> graph_ids = {std::numeric_limits<INDEX>::max(), std::numeric_limits<INDEX>::max()};
    std::pair<Nodes, Nodes> node_mapping = {std::vector<NodeId>(), std::vector<NodeId>()};
    std::string graph_data_name;
    double time = std::numeric_limits<double>::max();

    void get_edit_operations(std::unordered_set<EditOperation, EditOperationHash>& edit_operations) const;
    void get_edit_path(EditPath<T>& edit_path, int seed = 0) const;

    // edit path extension
    void delete_edge(EditPath<T>& edit_path, const EditOperation& operation) const;
    static void insert_edge(EditPath<T>& edit_path, const EditOperation& operation);
    static void relabel_edge(EditPath<T>& edit_path, const EditOperation& operation);
    void delete_node(EditPath<T>& edit_path, const EditOperation& operation) const;
    static void insert_node(EditPath<T>& edit_path, const EditOperation& operation);
    static void relabel_node(EditPath<T>& edit_path, const EditOperation& operation);

};

template <typename T>
inline void GEDEvaluation<T>::get_edit_operations(std::unordered_set<EditOperation, EditOperationHash>& edit_operations) const {
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
    for (auto s_id : mapping_second) {
        if (s_id > mapping_first.size()) {
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
                    if (!graphs.second.IsEdge(target_i, target_j)) {
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
                    if (!graphs.first.IsEdge(source_i, source_j)) {
                        inserted_edges.insert({target_i, target_j});
                    }
                } else {
                    inserted_edges.insert({target_i, target_j});
                }
            }
        }
    }
    // Find relabeled edges
    for (auto source_i = 0; source_i < graphs.first.nodes(); ++source_i) {
        for (auto source_j : graphs.first.get_neighbors(source_i)) {
            if (source_i < source_j) {
                NodeId target_i = mapping_first[source_i];
                NodeId target_j = mapping_first[source_j];
                if (target_i < mapping_second.size() && target_j < mapping_second.size()) {
                    if (graphs.first.IsEdge(source_i, source_j) && graphs.second.IsEdge(target_i, target_j)) {
                        // check if labels are different
                        if (graphs.first.GetEdgeData({source_i, source_j}, "label") != graphs.second.GetEdgeData({target_i, target_j}, "label")) {
                            relabeled_edges.insert({source_i, source_j});
                        }
                    }
                }
            }
        }
    }
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

template <typename T>
inline void GEDEvaluation<T>::get_edit_path(EditPath<T>& edit_path, int seed) const {
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

template <typename T>
inline void GEDEvaluation<T>::delete_edge(EditPath<T> &edit_path, const EditOperation &operation) const {
    const NodeId source_i = operation.edge.first;
    const NodeId source_j = operation.edge.second;
    const NodeId current_i = edit_path.source_to_current[source_i];
    const NodeId current_j = edit_path.source_to_current[source_j];
    // if an edge deletion creates an isolated node, delete the node as well
    // new graph after deleting the edge
    T new_graph = edit_path.edit_path_graphs.back();
    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);
    new_graph.RemoveEdge(current_i, current_j);
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

template <typename T>
inline void GEDEvaluation<T>::delete_node(EditPath<T> &edit_path, const EditOperation &operation) const {
    const NodeId source_node = operation.node;
    const NodeId current_node = edit_path.source_to_current[source_node];



    // TODO check whether node has degree 0 o.w. delete all edges
    const T& last_graph = edit_path.edit_path_graphs.back();
    if (last_graph.degree(current_node) != 0) {
        // raise an error
        std::cerr << "Error: Trying to delete a node that is not isolated. Node: " << source_node << " Current Target Node: " << current_node << std::endl;

    }
    // new graph
    T new_graph = last_graph;
    new_graph.RemoveNode(current_node);
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

    edit_path.edit_path_graphs.emplace_back(new_graph);
    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_node_deletions.erase(operation);
    edit_path.remaining_operations.erase(operation);
}

template <typename T>
inline void GEDEvaluation<T>::insert_edge(EditPath<T> &edit_path, const EditOperation &operation) {
    // Check if both nodes of the edge are already inserted
    const NodeId target_i = operation.edge.first;
    const NodeId target_j = operation.edge.second;

    T& last_graph = edit_path.edit_path_graphs.back();
    NodeId current_i = edit_path.target_to_current[target_i];
    NodeId current_j = edit_path.target_to_current[target_j];
    const bool first_node_inserted = current_i <= last_graph.nodes();
    const bool second_node_inserted = current_j <= last_graph.nodes();

    if (!first_node_inserted) {
        const EditOperation first_node_insertion = {
            .operationObject = OperationObject::NODE,
            .type = EditType::INSERT,
            .node = operation.edge.first,
        };
        // Update last graph if node is inserted
        last_graph = edit_path.edit_path_graphs.back();
        insert_node(edit_path, first_node_insertion);
    }
    if (!second_node_inserted) {
        const EditOperation second_node_insertion = {
            .operationObject = OperationObject::NODE,
            .type = EditType::INSERT,
            .node = operation.edge.second,
        };
        insert_node(edit_path, second_node_insertion);
        // Update last graph if node is inserted
        last_graph = edit_path.edit_path_graphs.back();
    }

    // new graph after inserting the edge
    current_i = edit_path.target_to_current[target_i];
    current_j = edit_path.target_to_current[target_j];
    T new_graph = last_graph;
    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);
    // TODO edge features
    std::vector<double> edge_features = std::vector<double>(0);
    new_graph.AddEdge(current_i, current_j, edge_features, true);
    edit_path.edit_path_graphs.emplace_back(new_graph);

    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_operations.erase(operation);
    edit_path.remaining_edge_insertions.erase(operation);
}

template <typename T>
inline void GEDEvaluation<T>::insert_node(EditPath<T> &edit_path, const EditOperation &operation) {
    const T& last_graph = edit_path.edit_path_graphs.back();
    // new graph after inserting the node
    auto new_graph = last_graph;
    if (new_graph.labelType != LABEL_TYPE::UNLABELED) {
        new_graph.AddNodes(1, {edit_path.target_graph.label(operation.node)});
    } else {
        new_graph.AddNodes(1);
    }
    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);


    //Update the maps
    edit_path.target_to_current[operation.node] = last_graph.nodes();

    edit_path.edit_path_graphs.emplace_back(new_graph);
    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_node_insertions.erase(operation);
    edit_path.remaining_operations.erase(operation);
}

template <typename T>
inline void GEDEvaluation<T>::relabel_edge(EditPath<T> &edit_path, const EditOperation &operation) {
    const T& last_graph = edit_path.edit_path_graphs.back();
    T new_graph = last_graph;
    std::vector<double> target_edge_data = edit_path.target_graph.GetEdgeData(EDGE(operation.edge.first, operation.edge.second));
    new_graph.RelabelEdge(operation.edge.first, operation.edge.second, target_edge_data);

    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);

    edit_path.edit_path_graphs.emplace_back(new_graph);
    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_operations.erase(operation);
    edit_path.remaining_node_relabels.erase(operation);
}

template <typename T>
inline void GEDEvaluation<T>::relabel_node(EditPath<T> &edit_path, const EditOperation &operation) {
    const NodeId source_i = operation.node;
    const NodeId current_i = edit_path.source_to_current[source_i];
    const NodeId target_i = edit_path.node_mapping.first[source_i];

    const T& last_graph = edit_path.edit_path_graphs.back();
    T new_graph = last_graph;
    Label target_label = edit_path.target_graph.label(target_i);
    new_graph.RelabelNode(current_i, target_label);

    const std::string name = edit_path.source_graph.GetName() + "_step_" + std::to_string(edit_path.edit_path_graphs.size()) + "_" + edit_path.target_graph.GetName();
    new_graph.SetName(name);

    edit_path.edit_path_graphs.emplace_back(new_graph);
    edit_path.sequence_of_operations.push_back(operation);
    edit_path.remaining_operations.erase(operation);
    edit_path.remaining_node_relabels.erase(operation);
}

#endif //GEDEXAMPLE_GEDEVALUATION_H