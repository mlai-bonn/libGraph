//
// Created by florian on 15.12.21.
//

#ifndef HOPS_SIMPLEPATTERNS_H
#define HOPS_SIMPLEPATTERNS_H


#include <list>
#include "DataClasses.h"

class SimplePatterns {
public:
    static GraphStruct Circle(int size, const Labels* labels = nullptr, int num_labels = 2);
    static GraphStruct StarGraph(int size, const Labels* labels = nullptr, int num_labels = 2);
    static GraphStruct Triangle(const Labels* labels = nullptr, int num_labels = 2);

    static GraphStruct FullyConnected(int i, const Labels* labels = nullptr, int num_labels = 2);

    static GraphStruct DoubleTriangle(const Labels *labels = nullptr, int num_labels = 2);

    static GraphStruct Path(int i, Labels *labels = nullptr, int num_labels = 2);
    static GraphStruct ErdosRenyi(int size, int edges, int seed = 0, bool connected = false);
};

inline GraphStruct SimplePatterns::Triangle(const Labels* labels, int num_labels) {
    return Circle(3, labels, num_labels);
}

inline GraphStruct SimplePatterns::DoubleTriangle(const Labels* labels, int num_labels) {
    GraphStruct G = GraphStruct(4, Labels());
    G.add_edge(0, 1);
    G.add_edge(0, 2);
    G.add_edge(1, 2);
    G.add_edge(1, 3);
    G.add_edge(2,3);
    G.SetName("double_triangle");
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.set_labels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }


    return G;
}

inline GraphStruct SimplePatterns::FullyConnected(int i, const Labels* labels, int num_labels) {
    GraphStruct G = GraphStruct(i, Labels());
    for (int j = 0; j < i; ++j) {
        for (int k = j + 1; k < i; ++k) {
            G.add_edge(j, k);
        }
    }
    G.SetName("fully_connected");
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.set_labels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() +std::to_string(l));
        }
    }
    return G;
}

inline GraphStruct SimplePatterns::Path(int i, Labels *labels, int num_labels) {
    GraphStruct G = GraphStruct(i + 1, Labels());
    for (int j = 0; j < i; ++j) {
        G.add_edge(j, j + 1);
    }
    G.SetName("path_" + std::to_string(i));

    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.set_labels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }
    return G;
}

inline GraphStruct SimplePatterns::Circle(int size, const Labels *labels, int num_labels) {
    GraphStruct G = GraphStruct(size, Labels());
    for (int i = 0; i < size; ++i) {
        G.add_edge(i, (i + 1) % size);
    }
    G.SetName("circle_" + std::to_string(size));
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.set_labels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }
    return G;
}

inline GraphStruct SimplePatterns::StarGraph(int size, const Labels *labels, int num_labels) {
    GraphStruct G = GraphStruct(size, Labels());
    for (int i = 0; i < size - 1; ++i) {
        G.add_edge(0, i + 1);
    }
    G.SetName("star_" + std::to_string(size));
    if (labels != nullptr){
        G.SetNumLabels(num_labels);
        G.set_labels(labels);
        G.SetName(G.GetName() + "_labels_");
        for (auto l : *labels) {
            G.SetName(G.GetName() + std::to_string(l));
        }
    }
    return G;
}

inline GraphStruct SimplePatterns::ErdosRenyi(int size, int edges, int seed, bool connected) {
    std::mt19937_64 gen(seed);
    GraphStruct G = GraphStruct(size, Labels());
    if (connected){
        int max_tries = 1000000;
        int tries = 0;
        while (!GraphStruct::IsConnected(G) && tries < max_tries) {
            G = GraphStruct(size, Labels());
            while (G.edges() < edges) {
                int src = std::uniform_int_distribution<int>(0, size - 1)(gen);
                int dst = std::uniform_int_distribution<int>(0, size - 1)(gen);
                if (src != dst) {
                    G.add_edge(src, dst);
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
                G.add_edge(src, dst);
            }
        }
    }
    return G;
}


#endif //HOPS_SIMPLEPATTERNS_H
