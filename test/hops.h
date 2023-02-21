//
// Created by florian on 23.01.23.
//

#ifndef TESTGRAPHLIB_HOPS_H
#define TESTGRAPHLIB_HOPS_H

#include "../include/libGraph.h"

void HopsTest(){
    GraphStruct triangle = SimplePatterns::Triangle();
    GraphStruct star = SimplePatterns::StarGraph(3);
    GraphData<GraphStruct> graph_data = GraphData<GraphStruct>();
    graph_data.add(SimplePatterns::DoubleTriangle());
    Hops hops = Hops(graph_data);
    RunParameters rParameters = {LABEL_TYPE::UNLABELED,0,1, 10000};
    hops.Run(0, star, rParameters);
}

void HopsAutomorphismTest(){
    for (int size = 2; size < 10; ++size) {
        GraphData<GraphStruct> graphData = GraphData<GraphStruct>("../../GraphData/Hops/patterns/size" + std::to_string(size) + ".bgfs");
        Hops hops = Hops(graphData);
        FileEvaluation fileEvaluation = FileEvaluation("../test/Results/", "automorphisms_" + std::to_string(size));
        for (int i = 0; i < graphData.size(); ++i) {
            hops.Automorphisms(i, {LABEL_TYPE::UNLABELED,30});
            fileEvaluation.headerValueInsert({"Size", "Edges", "Automorphisms"},{std::to_string(graphData[i].nodes()), std::to_string(graphData[i].edges()),
                                                                                 std::to_string(std::round(hops.hopsEvaluation.hopsEstimation))});
        }
        fileEvaluation.save();
    }
}

void HopsTimeTest(){
    GraphStruct triangle = SimplePatterns::Triangle();
    GraphStruct five_clique = SimplePatterns::FullyConnected(5);
    Hops hops = Hops("../../../../GraphData/Hops/com-amazon.ungraph.bgfs", "../Results/time_test_");
    hops.Run(0, triangle, {LABEL_TYPE::UNLABELED,0,1, 100000000, 0});
    hops.Run(0, five_clique, {LABEL_TYPE::UNLABELED,0,1, 100000000, 0});
}

void Hops60s(){
    GraphStruct triangle = SimplePatterns::Triangle();
    GraphData graphs = GraphData<GraphStruct>("../../../../GraphData/Hops/com-amazon.ungraph.bgfs");
    Hops hops = Hops(graphs);
    hops.Run(0, triangle, {LABEL_TYPE::UNLABELED,60,-1, 0, 0,true});
}

void HopsPatternTest(){
    std::cout << " Counting patterns for " << "amazon" << std::endl;
    for (auto size : {2,3,4,5,6}) {
        std::vector<std::vector<std::string>> out;
        StaticFunctionsLib::load_csv("../../GraphData/Hops/patterns/automorphisms_" + std::to_string(size) + ".csv", out);
        std::vector<int> automorphisms;
        int counter = 0;
        for (auto &x : out) {
            if (counter > 0) {
                automorphisms.emplace_back(std::stoi(x[3]));
            }
            ++counter;
        }

        Hops hops = Hops("../../GraphData/Hops/com-amazon.ungraph.bgfs", "../test/Results/pattern_test_");

        std::cout << std::endl << "Counting patterns of size " << size << std::endl;
        GraphData g_patterns = GraphData<GraphStruct>("../../GraphData/Hops/patterns/size" + std::to_string(size)+ ".bgfs");
        for (auto& pattern : g_patterns.graphData) {
            RunParameters runParameters = {LABEL_TYPE::UNLABELED, 5, -1, 0, 0, true, true, true, 1};
            hops.Run(0, pattern, runParameters);
        }
    }

}

void HopsParallelizationTest(){
    GraphStruct triangle = SimplePatterns::Triangle();

    for (auto graph : {"dblp", "amazon", "youtube", "orkut", "lj"}) {
        Hops hops = Hops("../../../../GraphData/Hops/com-" + graph + ".ungraph.bgfs", "../test/Results/");
        int max_threads = omp_get_max_threads();
        for (int i = 1; i <= max_threads; ++i) {
            RunParameters unlabeledRun{LABEL_TYPE::UNLABELED, 30, i, 0, 0, true, true, true};
            hops.Run(0, triangle, unlabeledRun);
        }
        for (int i = 1; i <= max_threads; ++i) {
            RunParameters unlabeledRun{LABEL_TYPE::UNLABELED, 0, i, 100000000/i, 0, true, true, true};
            hops.Run(0, triangle, unlabeledRun);
        }
    }
}

void HopsRealWorldTest() {
    GraphStruct triangle = SimplePatterns::Triangle();
    GraphStruct graph1 = GraphStruct("../../GraphData/Hops/com-amazon.ungraph.bgfs");
    GraphStruct graph2 = GraphStruct("../../GraphData/Hops/com-youtube.ungraph.bgfs");
    GraphData graph_data = GraphData<GraphStruct>();
    graph_data.add({graph1, graph2});
    RunParameters unlabeledRun{LABEL_TYPE::UNLABELED, 60, -1, 0, 0, true, true};


    Hops hops = Hops(graph_data);
    hops.Run(0, triangle, unlabeledRun);
    hops.Run(1, triangle, unlabeledRun);
}

void HopsSimplePatternsTest(){
    Labels triangleLabels = Labels{0, 0, 1};
    GraphStruct triangle = SimplePatterns::Triangle(&triangleLabels);

    Labels doubleTriangleLabels = Labels{1,0,1,0};
    GraphStruct doubleTriangle = SimplePatterns::DoubleTriangle(&doubleTriangleLabels);

    Labels fullyConnectedLabels = Labels(20, 0);
    fullyConnectedLabels[0] = 1;
    GraphStruct fullyConnected = SimplePatterns::FullyConnected(20, &fullyConnectedLabels);


    Labels path2Labels = Labels {1, 0, 1};
    GraphStruct path2 = SimplePatterns::Path(2, &path2Labels);

    Labels path3Labels = Labels{1, 0, 0, 1};
    GraphStruct path3 = SimplePatterns::Path(3, &path3Labels);


    Labels squareLabels1 = Labels{0, 1, 0, 1};
    Labels squareLabels2 = Labels{0, 1, 1, 0};
    GraphStruct square1 = SimplePatterns::Circle(4, &squareLabels1);
    GraphStruct square2 = SimplePatterns::Circle(4, &squareLabels2);

    GraphData graph_data = GraphData<GraphStruct>();
    graph_data.add({doubleTriangle, fullyConnected});

    RunParameters labeledRun{LABEL_TYPE::LABELED_SPARSE, 1, 1, 0, 0, true, false};
    RunParameters unlabeledRun{LABEL_TYPE::UNLABELED, 1, 1, 0, 0, true, false};


    Hops hops = Hops(graph_data);
    hops.Run(0, triangle, labeledRun);
    hops.Run(0, triangle, unlabeledRun);
    hops.Run(1, triangle, labeledRun);
    hops.Run(1, triangle, unlabeledRun);
    hops.Run(0, path2, labeledRun);
    hops.Run(0, path2, unlabeledRun);
    hops.Run(1, path2, labeledRun);
    hops.Run(1, path2, unlabeledRun);
    hops.Run(0, path3, labeledRun);
    hops.Run(0, path3, unlabeledRun);
    hops.Run(1, path3, labeledRun);
    hops.Run(1, path3, unlabeledRun);
    hops.Run(0, square1, labeledRun);
    hops.Run(0, square1, unlabeledRun);
    hops.Run(1, square1, labeledRun);
    hops.Run(1, square1, unlabeledRun);
    hops.Run(0, square2, labeledRun);
    hops.Run(0, square2, unlabeledRun);
    hops.Run(1, square2, labeledRun);
    hops.Run(1, square2, unlabeledRun);
}

#endif //TESTGRAPHLIB_HOPS_H
