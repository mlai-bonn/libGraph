//
// Created by florian on 17.10.23.
//

#ifndef LIBGRAPH_GRAPHALGORITHMS_H
#define LIBGRAPH_GRAPHALGORITHMS_H

#include <stack>

/**
 * @brief Computes the bi-connected components of a graph
 * @param graph The input graph
 * @param components The output components (each component is a vector of node ids)
 */
static void GetBiconnectedComponents(const GraphStruct& graph, std::vector<std::vector<NodeId>>& components);
static void GetBiconnectedComponents(const GraphStruct& graph, std::vector<GraphStruct>& components);

static void GetBiconnectedOuterplanarFaces(const GraphStruct &component, OuterplanarComponent& outerplanarComponent);
static void GetBiconnectedOuterplanarFaceNum(GraphStruct &component, int &face_num);



inline void GetBiconnectedComponents(const GraphStruct& graph, std::vector<std::vector<NodeId>>& components){
    // Algorithm of Hopcroft and Tarjan
    // https://en.wikipedia.org/wiki/Biconnected_component

    std::vector<NodeId> lowPoint(graph.nodes(), 0);
    std::vector<NodeId> depth(graph.nodes(), 0);
    std::stack<NodeId> stack;
    std::vector<bool> visited(graph.nodes(), false);
    std::vector<bool> backtracked(graph.nodes(), false);
    std::vector<bool> articulationPoints(graph.nodes(), false);
    std::vector<NodeId> leaves;

    std::vector<NodeId> parent(graph.nodes(), 0);
    components.clear();

    // initialize the root vertex
    stack.push(0);
    depth[0] = 0;
    parent[0] = -1;

    while (!stack.empty()){
        NodeId current_node = stack.top();
        if (!visited[current_node]) {
            visited[current_node] = true;
            int count_visited_neighbors = 0;
            for (auto neighbor : graph.get_neighbors(current_node)) {
                if (!visited[neighbor]) {
                    depth[neighbor] = depth[current_node] + 1;
                    lowPoint[neighbor] = depth[neighbor];
                    stack.push(neighbor);
                    parent[neighbor] = current_node;
                    ++count_visited_neighbors;
                }
            }
            if (count_visited_neighbors == 0) {
                leaves.emplace_back(current_node);
            }
        }
        else {
            if (!backtracked[current_node]) {
                backtracked[current_node] = true;
                // we are backtracking
                if (current_node != 0) {
                    // check the neighbors with lower depth
                    for (auto neighbor: graph.get_neighbors(current_node)) {
                        if (neighbor != parent[current_node]) {
                            lowPoint[current_node] = std::min(lowPoint[current_node], depth[neighbor]);
                            if (depth[neighbor] - 1 == depth[current_node]) {
                                lowPoint[current_node] = std::min(lowPoint[current_node], lowPoint[neighbor]);
                            }
                            if (lowPoint[neighbor] >= depth[current_node]) {
                                articulationPoints[current_node] = true;
                                break;
                            }
                        }
                    }
                }
                // handle the root vertex separately
                else {
                    int children = 0;
                    for (auto neighbor: graph.get_neighbors(current_node)) {
                        if (depth[neighbor] == 1) {
                            ++children;
                        }
                        if (children > 1) {
                            articulationPoints[current_node] = true;
                            break;
                        }
                    }
                }
            }
            stack.pop();
        }
    }
    // compute the components using the leaves the parents and the articulation points
    std::fill(visited.begin(), visited.end(), false);
    components.emplace_back();
    while (!leaves.empty()) {
        NodeId current_node = leaves.back();
        leaves.pop_back();
        if (current_node != 0 && !visited[current_node]) {
            if (articulationPoints[current_node] && articulationPoints[parent[current_node]]){
                components.back().emplace_back(current_node);
                visited[current_node] = true;
                components.back().emplace_back(parent[current_node]);
                visited[current_node] = true;
                components.emplace_back();
                current_node = parent[current_node];
                if (!visited[current_node]) {
                    leaves.emplace_back(current_node);
                }
            }
            else if (articulationPoints[current_node] && !articulationPoints[parent[current_node]]){
                components.back().emplace_back(current_node);
                visited[current_node] = true;
                current_node = parent[current_node];
                if (!visited[current_node]) {
                    leaves.emplace_back(current_node);
                }
            }
            else{
                while (current_node != 0 && !articulationPoints[current_node]) {
                    visited[current_node] = true;
                    components.back().emplace_back(current_node);
                    current_node = parent[current_node];
                    if (articulationPoints[current_node]){
                        if (!visited[current_node]) {
                            leaves.emplace_back(current_node);

                            leaves.emplace_back(current_node);
                            components.back().emplace_back(current_node);
                            components.emplace_back();
                        }
                        else{
                            components.back().emplace_back(current_node);
                        }
                    }
                }
            }
        }
    }
    return;
}

void GetBiconnectedComponents(const GraphStruct &graph, std::vector<GraphStruct> &components) {
    std::vector<std::vector<NodeId>> component_nodes;
    GetBiconnectedComponents(graph, component_nodes);
    components.clear();
    for (const auto& component : component_nodes){
        if (component.size() > 2) {
            components.emplace_back(GraphStruct::SubGraph(graph, component));
        }
    }
}

inline void GetBiconnectedOuterplanarFaces(const GraphStruct &component, OuterplanarComponent &outerplanarComponent) {
    std::vector<std::vector<NodeId>> neighbors = std::vector<std::vector<NodeId>>(component.nodes());
    std::vector<GraphStruct> currentFace;
    int face_num = 0;
    Nodes degree2Nodes;
    std::vector<int> degrees = std::vector<int>(component.nodes(), 0);
    //Preprocessing on get_node degrees
    for (NodeId id = 0; id < component.nodes(); ++id) {
        INDEX degree = component.degree(id);
        if (degree == 2) {
            degree2Nodes.emplace_back(id);
        }
        neighbors[id] = {component.neighbor(id, 0), component.neighbor(id, 1)};
        degrees[id] = degree;
    }
    outerplanarComponent.faces.emplace_back();
    currentFace.emplace_back(outerplanarComponent.faces.back());
    while (!degree2Nodes.empty() && component.edges() > 2){
        NodeId currentNodeId = degree2Nodes.back();
        degree2Nodes.pop_back();
        if (degrees[currentNodeId] > 0) {
            degrees[currentNodeId] = 0;
            NodeId n1 = neighbors[currentNodeId][0];
            NodeId n2 = neighbors[currentNodeId][1];

            int near_degree = degrees[n1];
            int next_degree = degrees[n2];
            if (component.edge(n1, n2)) {
                NodeId newNode1 = currentFace.back().add_node();
                NodeId newNode2 = currentFace.back().add_node();
                currentFace.back().add_edge_no_check(newNode1, newNode2);
                ++face_num;
                currentFace.pop_back();
                //Update degreeTwo nodes
                --degrees[n1];
                --degrees[n2];
            } else {
                NodeId newCurrent = currentFace.back().add_node();
                NodeId newNode1 = currentFace.back().add_node();
                NodeId newNode2 = currentFace.back().add_node();

                currentFace.back().add_edge_no_check(newCurrent, newNode1);
                currentFace.back().add_edge_no_check(newCurrent, newNode2);
            }
            if (degrees[n1] == 2 || degrees[n2] == 2) {
                if (degrees[n1] == 2) {
                    degree2Nodes.emplace_back(n1);
                    for (int i = 0; i < component.degree(n1); ++i) {
                        NodeId neighbor = component.neighbor(n1, i);
                        if (i > 0) {
                            neighbors[n1].emplace_back(neighbor);
                        }
                        if (neighbors[n1].size() == 2) {
                            break;
                        }
                    }
                } else {
                    degree2Nodes.emplace_back(n2);
                    for (int i = 0; i < component.degree(n2); ++i) {
                        NodeId neighbor = component.neighbor(n2, i);
                        if (i > 0) {
                            neighbors[n2].emplace_back(neighbor);
                        }
                        if (neighbors[n2].size() == 2) {
                            break;
                        }
                    }
                }
            }
            else {
                outerplanarComponent.faces.emplace_back();
                currentFace.emplace_back(outerplanarComponent.faces.back());
            }
        }
    }
}

inline void GetBiconnectedOuterplanarFaceNum(GraphStruct &component, int &face_num) {
    face_num = 0;
    Nodes degree2Nodes;
    //Preprocessing on get_node degrees
    for (NodeId id = 0; id < component.nodes(); ++id) {
        INDEX degree = component.degree(id);
        if (degree == 2) {
            degree2Nodes.emplace_back(id);
        }
    }
    while (!degree2Nodes.empty() && component.edges() > 2){
        NodeId currentNodeId = degree2Nodes.back();
        degree2Nodes.pop_back();
        NodeId n1 = component.neighbor(currentNodeId,0);
        NodeId n2 = component.neighbor(currentNodeId,1);
        //Delete possibleEdges from component
        component.remove_edge(currentNodeId, n1);
        component.remove_edge(currentNodeId, n2);
        if (component.edge(n1, n2)) {
            ++face_num;
            //Update degreeTwo nodes
            INDEX deg1 = component.degree(n1);
            INDEX deg2 = component.degree(n2);
            if (deg1 == 2){
                degree2Nodes.emplace_back(n1);
            }
            if (deg2 == 2){
                degree2Nodes.emplace_back(n2);
            }
        }
        else{
            component.add_edge(n1, n2);
        }
    }
}


#endif //LIBGRAPH_GRAPHALGORITHMS_H
