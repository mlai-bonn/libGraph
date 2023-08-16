//
// Created by florian on 16.08.23.
//

#ifndef TESTGRAPHLIB_CLOSURES_H
#define TESTGRAPHLIB_CLOSURES_H

#include "../include/DataClasses.h"
#include "../include/Closures/GraphClosures.h"
#include "../include/io/StaticFunctions.h"

void TestClosure(){
    // load a graph
    GraphStruct graph = GraphStruct("../../../../GraphData/RealWorld/Collaboration/CA-GrQc_component.edges");
    // create a closure object
    GraphClosureSP closure = GraphClosureSP(graph);
    // define the closure parameters
    ClosureParameters parameters = {.input_set =  {0, 1, 5, 6}};
    // compute the closure
    closure.closure(parameters);
    // print the result
    std::cout << "Closure: " << StaticFunctionsLib::print<std::set<NodeId>, NodeId>(parameters.closed_set) << std::endl;

}

#endif //TESTGRAPHLIB_CLOSURES_H
