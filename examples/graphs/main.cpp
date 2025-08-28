// This file shows how to use the libGraph graph library


#include <iostream>
#include <string>
#include "GraphDataStructures/GraphBase.h"

// Creating Simple (Unlabeled) Graphs
void unlabeled_graphs() {

    // Triangle
    std::string name = "Triangle";
    auto triangle = GraphStruct(name);
    // Adding nodes
    for (int i = 0; i < 3; ++i) {
        triangle.add_node();
    }
    // Adding edges
    triangle.add_edge(0, 1);
    triangle.add_edge(1, 2);
    triangle.add_edge(2, 0);
    // Printing the graph
    std::cout << triangle << std::endl;
}

void labeled_graphs() {

    //Triangle
    std::string name = "Triangle";
    auto triangle = GraphStruct(name);
    // Adding nodes
    for (int i = 0; i < 3; ++i) {
        triangle.add_node();
    }
    const Labels labels{1,1,2};
    triangle.set_labels(&labels);
    // Adding edges
    triangle.add_edge(0, 1);
    triangle.add_edge(1, 2);
    triangle.add_edge(2, 0);
    // Printing the graph
    std::cout << triangle << std::endl;

}

// main function to run the examples

int main() {
    unlabeled_graphs();
    labeled_graphs();
    return 0;

}

