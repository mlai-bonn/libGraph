//
// Created by florian on 17.10.23.
//

#ifndef LIBGRAPH_GRAPHALGORITHMS_H
#define LIBGRAPH_GRAPHALGORITHMS_H

#include <stack>
#include "GraphDataStructures/GraphStructs.h"

/**
 * @brief Computes the bi-connected components of a _graph
 * @param graph The input _graph
 * @param components The output components (each component is a vector of node ids)
 */
static void GetBiconnectedComponents(const GraphStruct& graph, std::vector<std::vector<NodeId>>& components);
static void GetBiconnectedComponents(const GraphStruct& graph, std::vector<GraphStruct>& components);

static void GetBiconnectedOuterplanarFaces(const GraphStruct &component, OuterplanarComponent& outerplanarComponent);
static void GetBiconnectedOuterplanarFaceNum(GraphStruct &component, int &face_num);

static void CheckingOuterpanarity(const GraphStruct &graph, const GraphStruct &outerplanarSubgraph, int &notOuterplanarSubgraphs,
                                  int &notMaximalSubgraphs, std::vector<int> &nonOuterplanarSeeds,
                                  std::vector<int> &nonMaximalSeeds, std::vector<double> &algorithmMissingEdges,
                                  std::vector<double> &maximalEdges, int seed = 0);

static bool IsOuterPlanar(const GraphStruct &graph, NodeId src=-1, NodeId dst=-1);
static bool IsMaximalOuterplanarSubgraph(const GraphStruct &graph, const GraphStruct &subgraph, std::vector<NodePair> &missingEdges);


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
                currentFace.back().add_edge(newNode1, newNode2, false);
                ++face_num;
                currentFace.pop_back();
                //Update degreeTwo nodes
                --degrees[n1];
                --degrees[n2];
            } else {
                NodeId newCurrent = currentFace.back().add_node();
                NodeId newNode1 = currentFace.back().add_node();
                NodeId newNode2 = currentFace.back().add_node();

                currentFace.back().add_edge(newCurrent, newNode1, false);
                currentFace.back().add_edge(newCurrent, newNode2, false);
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

void CheckingOuterpanarity(const GraphStruct &graph, const GraphStruct &outerplanarSubgraph, int &notOuterplanarSubgraphs,
                           int &notMaximalSubgraphs, std::vector<int> &nonOuterplanarSeeds,
                           std::vector<int> &nonMaximalSeeds, std::vector<double> &algorithmMissingEdges,
                           std::vector<double> &maximalEdges, int seed) {
//GraphFunctions::print(_subgraph);
    bool outerplanar = IsOuterPlanar(outerplanarSubgraph);
    std::vector<NodePair> missingEdges;
    bool maximalOuterplanar = IsMaximalOuterplanarSubgraph(graph, outerplanarSubgraph, missingEdges);
    std::string answer;
    std::string maximalAnswer;
    if (outerplanar) {
        answer = "Yes";
    } else {
        ++notOuterplanarSubgraphs;
        nonOuterplanarSeeds.emplace_back(seed);
        answer = "No";
    }
    if (maximalOuterplanar) {
        maximalAnswer = "Yes";

    } else {
        ++notMaximalSubgraphs;
        nonMaximalSeeds.emplace_back(seed);
        maximalAnswer = "No";
    }
    algorithmMissingEdges.emplace_back(missingEdges.size());
//    if (printLevel == PrintLevel::ALL) {
//        std::cout << "Outerplanar: " << answer << ", maximal: " << maximalAnswer << std::endl;
//
//        if (!maximalOuterplanar) {
//            std::cout << "Missing Edges ";
//            for (auto const &edge: missingEdges) {
//                std::cout << "(" << edge.first() << ", " << edge.second() << ")";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "//////////////////////////////////////" << std::endl;
//    }
    maximalEdges.emplace_back(outerplanarSubgraph.edges() + missingEdges.size());
}

bool IsOuterPlanar(const GraphStruct &graph, NodeId src, NodeId dst) {
    std::vector<std::vector<NodeId>> components;
    GetBiconnectedComponents(graph, components);
    for (const auto& component : components) {
        if (component.size() > 2) {
            if (src == -1 || dst == -1 || (std::find(component.begin(), component.end(), src) != component.end() && std::find(component.begin(), component.end(), src) != component.end())) {
                GraphStruct componentGraph = GraphStruct::SubGraph(graph, component);
                if (componentGraph.edges() > 2 * componentGraph.nodes() - 3) {
                    return false;
                } else {
                    GraphStruct algorithmGraph;
                    Nodes degree2Nodes;
                    std::unordered_map<NodePair, int, hashNodePair> triangulationCount;
                    std::vector<NodePair> pairs;
                    std::vector<NodePair> edges;
                    //Preprocessing on get_node degrees
                    for (auto node : componentGraph) {
                        algorithmGraph.add_node();
                        INDEX degree = componentGraph.degree(node);
                        if (degree == 2) {
                            degree2Nodes.emplace_back(node);
                        }
                    }
                    for (auto edge = componentGraph.first_edge(); edge != componentGraph.last_edge(); ++edge) {
                        algorithmGraph.add_edge(*edge);
                        triangulationCount[NodePair(*edge, false)] = 0;
                    }
                    while (!degree2Nodes.empty()) {
                        NodeId cNode = degree2Nodes.back();
                        degree2Nodes.pop_back();
                        if (algorithmGraph.degree(cNode) == 2) {
                            NodeId near = algorithmGraph[cNode][0];
                            NodeId next = algorithmGraph[cNode][1];
                            NodePair nodePair = NodePair(near, next);
                            algorithmGraph.remove_edge(cNode, near);
                            algorithmGraph.remove_edge(cNode, next);
                            if (!algorithmGraph.edge(near, next)) {
                                algorithmGraph.add_edge(near, next);
                                triangulationCount[nodePair] = std::max(1, std::max(
                                        triangulationCount[NodePair(cNode, near)],
                                        triangulationCount[NodePair(cNode, next)]));
                                if (triangulationCount[nodePair] > 2) {
                                    return false;
                                }
                            } else {
                                triangulationCount[nodePair] += 1;
                                if (triangulationCount[nodePair] > 2 ||
                                    triangulationCount[NodePair(cNode, near)] == 2 ||
                                    triangulationCount[NodePair(cNode, next)] == 2) {
                                    return false;
                                }
                                if (algorithmGraph.degree(near) == 2) {
                                    degree2Nodes.emplace_back(near);
                                }
                                if (algorithmGraph.degree(next) == 2) {
                                    degree2Nodes.emplace_back(next);
                                }
                            }
                        }
                    }
                    if (algorithmGraph.edges() > 1) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool IsMaximalOuterplanarSubgraph(const GraphStruct &graph, const GraphStruct & subgraph, std::vector<NodePair>& missingEdges) {
    bool outerplanar = IsOuterPlanar(subgraph);
    if (!outerplanar){
        return false;
    }
    bool additionalEdge = false;
    missingEdges.clear();

    GraphStruct check_graph = subgraph;

    for (auto edge = graph.first_edge(); edge != graph.last_edge(); ++edge) {
        if (!check_graph.edge(*edge)){
            check_graph.add_edge(*edge);
            additionalEdge = IsOuterPlanar(check_graph, (*edge).first, (*edge).second);
            if (additionalEdge){
                missingEdges.emplace_back(*edge, false);
            }
            else{
                check_graph.remove_edge(*edge);
            }
        }
    }
    if (missingEdges.empty()){
        return true;
    }
    return false;
}



#endif //LIBGRAPH_GRAPHALGORITHMS_H
