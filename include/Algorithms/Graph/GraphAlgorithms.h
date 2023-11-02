//
// Created by florian on 17.10.23.
//

#ifndef LIBGRAPH_GRAPHALGORITHMS_H
#define LIBGRAPH_GRAPHALGORITHMS_H

#include <stack>
#include "GraphDataStructures/GraphStructs.h"

class GraphAlgorithms {
public:

/**
 * @brief Computes the bi-connected components of a _graph
 * @param graph The input _graph
 * @param components The output components (each component is a vector of node ids)
 */
    static void GetBiconnectedComponents(const GraphStruct &graph, std::vector<std::vector<NodeId>> &components);

    static void GetBiconnectedComponents(const GraphStruct &graph, std::vector<GraphStruct> &components);

    static void
    GetBiconnectedOuterplanarFaces(const GraphStruct &component, OuterplanarComponent &outerplanarComponent);

    static void GetBiconnectedOuterplanarFaceNum(GraphStruct &component, int &face_num);

    static void CheckingOuterpanarity(const GraphStruct &graph, const GraphStruct &outerplanarSubgraph,
                                      int &notOuterplanarSubgraphs,
                                      int &notMaximalSubgraphs, std::vector<int> &nonOuterplanarSeeds,
                                      std::vector<int> &nonMaximalSeeds, std::vector<double> &algorithmMissingEdges,
                                      std::vector<double> &maximalEdges, int seed = 0);

    static bool IsOuterPlanar(const GraphStruct &graph, NodeId src = -1, NodeId dst = -1);

    static bool IsMaximalOuterplanarSubgraph(const GraphStruct &graph, const GraphStruct &subgraph,
                                             std::vector<NodePair> &missingEdges);

    static void BCCUtil(const GraphStruct &graph, std::vector<std::vector<NodeId>> &components, NodeId u, std::vector<NodeId> &disc,
                 std::vector<NodeId> &low, std::vector<std::pair<NodeId, NodeId>> &edges, std::vector<NodeId> &parents,
                 std::vector<NodeId> &visitedId, int &discovery_time);

    static void AddRandomEdges(GraphStruct &graph, INDEX edgeNum, int seed);
};

inline void GraphAlgorithms::GetBiconnectedComponents(const GraphStruct &graph, std::vector<GraphStruct> &components) {
    std::vector<std::vector<NodeId>> component_nodes;
    GraphAlgorithms::GetBiconnectedComponents(graph, component_nodes);
    components.clear();
    for (const auto& component : component_nodes){
            components.emplace_back(GraphStruct::SubGraph(graph, component));
    }
}

inline void GraphAlgorithms::GetBiconnectedOuterplanarFaces(const GraphStruct &component, OuterplanarComponent &outerplanarComponent) {
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

inline void GraphAlgorithms::GetBiconnectedOuterplanarFaceNum(GraphStruct &component, int &face_num) {
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

void GraphAlgorithms::CheckingOuterpanarity(const GraphStruct &graph, const GraphStruct &outerplanarSubgraph, int &notOuterplanarSubgraphs,
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

bool GraphAlgorithms::IsOuterPlanar(const GraphStruct &graph, NodeId src, NodeId dst) {
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

bool GraphAlgorithms::IsMaximalOuterplanarSubgraph(const GraphStruct &graph, const GraphStruct & subgraph, std::vector<NodePair>& missingEdges) {
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


// A C++ program to find biconnected components in a given undirected graph

// A recursive function that finds and prints strongly connected
// components using DFS traversal
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
// discovery time) that can be reached from subtree
// rooted with current vertex
// *st -- >> To store visited edges
void GraphAlgorithms::BCCUtil(const GraphStruct& graph, std::vector<std::vector<NodeId>>& components, NodeId u, std::vector<NodeId>& disc, std::vector<NodeId>& low, std::vector<std::pair<NodeId, NodeId>>& edges, std::vector<NodeId>& parents, std::vector<NodeId>& visitedId, int& discovery_time)
{

    // Initialize discovery time and low value
    disc[u] = low[u] = ++discovery_time;
    int children = 0;

    // Go through all vertices adjacent to this
    for (NodeId i = 0; i < graph.degree(u); ++i) {
        NodeId v = graph.neighbor(u, i); // v is current adjacent of 'u'

        // If v is not visited yet, then recur for it
        if (disc[v] == -1) {
            children++;
            parents[v] = u;
            // store the edge in stack
            edges.emplace_back(u, v);
            BCCUtil(graph, components, v, disc, low, edges, parents, visitedId, discovery_time);

            // Check if the subtree rooted with 'v' has a
            // connection to one of the ancestors of 'u'
            // Case 1 -- per Strongly Connected Components Article
            low[u] = std::min(low[u], low[v]);

            // If u is an articulation point,
            // pop all edges from stack till u -- v
            if ((disc[u] == 1 && children > 1) || (disc[u] > 1 && low[v] >= disc[u])) {
                components.emplace_back();
                while (edges.back().first != u || edges.back().second != v) {
                    if (visitedId[edges.back().first] != discovery_time) {
                        components.back().emplace_back(edges.back().first);
                        visitedId[edges.back().first] = discovery_time;
                    }
                    if (visitedId[edges.back().second] != discovery_time) {
                        components.back().emplace_back(edges.back().second);
                        visitedId[edges.back().second] = discovery_time;
                    }
                    edges.pop_back();
                }
                if (visitedId[edges.back().first] != discovery_time) {
                    components.back().emplace_back(edges.back().first);
                    visitedId[edges.back().first] = discovery_time;
                }
                if (visitedId[edges.back().second] != discovery_time) {
                    components.back().emplace_back(edges.back().second);
                    visitedId[edges.back().second] = discovery_time;
                }
                edges.pop_back();
                discovery_time++;
            }
        }

        // Update low value of 'u' only of 'v' is still in stack
        // (i.e. it's a back edge, not cross edge).
        // Case 2 -- per Strongly Connected Components Article
        else if (v != parents[u]) {
            low[u] = std::min(low[u], disc[v]);
            if (disc[v] < disc[u]) {
                edges.emplace_back(u, v);
            }
        }
    }
}

// The function to do DFS traversal. It uses BCCUtil()
inline void GraphAlgorithms::GetBiconnectedComponents(const GraphStruct& graph, std::vector<std::vector<NodeId>>& components)
{
    components.clear();
    std::vector<NodeId> disc = std::vector<NodeId>(graph.nodes(), -1);
    std::vector<NodeId> low = std::vector<NodeId>(graph.nodes(), -1);
    std::vector<NodeId> parents = std::vector<NodeId>(graph.nodes(), -1);
    std::vector<std::pair<NodeId, NodeId>> edges;
    std::vector<NodeId> visitedId = std::vector<NodeId>(graph.nodes(), -1);
    int discovery_time = 0;

    for (int i = 0; i < graph.nodes(); i++) {
        if (disc[i] == -1) {
            GraphAlgorithms::BCCUtil(graph, components, i, disc, low, edges, parents, visitedId, discovery_time);
        }
        if (!edges.empty()) {
            int j = 0;
            components.emplace_back();
            // If stack is not empty, pop all edges from stack
            while (!edges.empty()) {
                j = 1;
                if (visitedId[edges.back().first] != discovery_time) {
                    components.back().emplace_back(edges.back().first);
                    visitedId[edges.back().first] = discovery_time;
                }
                if (visitedId[edges.back().second] != discovery_time) {
                    components.back().emplace_back(edges.back().second);
                    visitedId[edges.back().second] = discovery_time;
                }
                edges.pop_back();
            }
            if (j == 1) {
                discovery_time++;
            }
        }
    }
}

inline void GraphAlgorithms::AddRandomEdges(GraphStruct &graph, INDEX edgeNum, int seed) {
    std::mt19937_64 gen(seed);
    if (edgeNum > (graph.nodes() * graph.nodes() - 1) / 2 - graph.edges()){
        throw std::range_error("Cannot add " + std::to_string(edgeNum) + " edges! Number of edges must be smaller than " + std::to_string((graph.nodes() * graph.nodes() - 1) / 2 - graph.edges()));
    }
    for (INDEX i = 1; i < edgeNum + 1; ++i) {
        NodeId src = std::uniform_int_distribution<NodeId>(0, graph.nodes() - 1)(gen);
        NodeId dst = std::uniform_int_distribution<NodeId>(0, graph.nodes() - 1)(gen);
        if (src != dst && graph.add_edge(src, dst)){
            continue;
        }
        else{
            --i;
        }
    }
}

#endif //LIBGRAPH_GRAPHALGORITHMS_H
