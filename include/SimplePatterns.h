//
// Created by florian on 15.12.21.
//

#ifndef SIMPLEPATTERNS_H
#define SIMPLEPATTERNS_H


#include <list>
#include "GraphDataStructures/GraphBase.h"

class SimplePatterns {
public:
    template<typename T>
    static T Circle(int size, const Labels* labels = nullptr, int num_labels = 2);
    template<typename T>
    static T StarGraph(int size, const Labels* labels = nullptr, int num_labels = 2);
    template<typename T>
    static T Triangle(const Labels* labels = nullptr, int num_labels = 2);
    template<typename T>
    static T MaxPreClosure(int i, Labels *labels = nullptr, int num_labels = 2);

    template<typename T>
    static T FullyConnected(int i, const Labels* labels = nullptr, int num_labels = 2);

    template<typename T>
    static T DoubleTriangle(const Labels *labels = nullptr, int num_labels = 2);

    template<typename T>
    static T Path(int i, Labels *labels = nullptr, int num_labels = 2);
    template<typename T>
    static T ErdosRenyi(int size, int edges, int seed = 0, bool connected = false);

    template<typename T>
    static T FullyBipartite(int partitionASize, int partitionBSize);
};

template<typename T>
inline T SimplePatterns::Triangle(const Labels* labels, int num_labels) {
    return Circle<T>(3, labels, num_labels);
}

template<typename T>
inline T SimplePatterns::DoubleTriangle(const Labels* labels, int num_labels) {
    T G = T(4, Labels());
    G.AddEdge(0, 1);
    G.AddEdge(0, 2);
    G.AddEdge(1, 2);
    G.AddEdge(1, 3);
    G.AddEdge(2,3);
    G.SetName("double_triangle");
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.SetLabels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }


    return G;
}

template<typename T>
inline T SimplePatterns::FullyConnected(int i, const Labels* labels, int num_labels) {
    T G = T(i, Labels());
    for (int j = 0; j < i; ++j) {
        for (int k = j + 1; k < i; ++k) {
            G.AddEdge(j, k);
        }
    }
    G.SetName("fully_connected");
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.SetLabels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() +std::to_string(l));
        }
    }
    return G;
}

template<typename T>
inline T SimplePatterns::Path(int i, Labels *labels, int num_labels) {
    T G = T(i + 1, Labels());
    for (int j = 0; j < i; ++j) {
        G.AddEdge(j, j + 1);
    }
    G.SetName("path_" + std::to_string(i));
    G.SetType(GraphType::TREE);

    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.SetLabels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }
    return G;
}

template<typename T>
inline T SimplePatterns::Circle(int size, const Labels *labels, int num_labels) {
    T G = T(size, Labels());
    for (int i = 0; i < size; ++i) {
        G.AddEdge(i, (i + 1) % size);
    }
    G.SetName("circle_" + std::to_string(size));
    G.SetType(GraphType::OUTERPLANAR);
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.SetLabels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }
    return G;
}

template<typename T>
inline T SimplePatterns::StarGraph(int size, const Labels *labels, int num_labels) {
    T G = T(size, Labels());
    for (int i = 0; i < size - 1; ++i) {
        G.AddEdge(0, i + 1);
    }
    G.SetName("star_" + std::to_string(size));
    G.SetType(GraphType::TREE);
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.SetLabels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }
    return G;
}

template<typename T>
inline T SimplePatterns::ErdosRenyi(int size, int edges, int seed, bool connected) {
    std::mt19937_64 gen(seed);
    T G = T(size, Labels());
    if (connected){
        int max_tries = 1000000;
        int tries = 0;
        while (!T::IsConnected(G) && tries < max_tries) {
            G = T(size, Labels());
            while (G.edges() < edges) {
                int src = std::uniform_int_distribution<int>(0, size - 1)(gen);
                int dst = std::uniform_int_distribution<int>(0, size - 1)(gen);
                if (src != dst) {
                    G.AddEdge(src, dst);
                }
            }
            ++tries;
        }
    }
    else {
        while (G.edges() < edges) {
            int src = std::uniform_int_distribution<int>(0, size - 1)(gen);
            int dst = std::uniform_int_distribution<int>(0, size - 1)(gen);
            if (src != dst) {
                G.AddEdge(src, dst);
            }
        }
    }
    G.SetName("erdos_renyi_" + std::to_string(size) + "_" + std::to_string(edges));
    return G;
}

template<typename T>
T SimplePatterns::MaxPreClosure(int i, Labels *labels, int num_labels) {
    T G = T(4+2*i, Labels());
    G.AddEdge(0, 2);
    G.AddEdge(0, 3);
    G.AddEdge(1, 2);
    G.AddEdge(1, 3);
    for (int j = 0; j < i; ++j) {
        G.AddEdge(2*j + 2, 2*j + 4);
        G.AddEdge(2*j + 2, 2*j + 5);
        G.AddEdge(2*j + 3, 2*j + 4);
        G.AddEdge(2*j + 3, 2*j + 5);
    }
    G.SetName("max_pre_closure_" + std::to_string(i));
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.SetLabels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }
    return G;
}

template<typename T>
T SimplePatterns::FullyBipartite(int partitionASize, int partitionBSize) {
    T G = T(partitionASize + partitionBSize, Labels());
    for (int i = 0; i < partitionASize; ++i) {
        for (int j = partitionASize; j < partitionASize + partitionBSize; ++j) {
            G.AddEdge(i, j);
        }
    }
    G.SetName("fully_bipartite_" + std::to_string(partitionASize) + "_" + std::to_string(partitionBSize));
    return G;
}


#endif //SIMPLEPATTERNS_H
