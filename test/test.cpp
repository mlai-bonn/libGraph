
#include "load_save.h"
#include "graph_functions.h"
#include "hops.h"

int main(int argc, char *argv[]) {

    int threads = -1;

    for (int i = 0; i < argc; ++i) {
        std::string str = argv[i];
        if (std::strcmp(argv[i], "--threads") == 0 || std::strcmp(argv[i], "--all")==0){
            threads = std::stoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--load_speed") == 0 || std::strcmp(argv[i], "--all") == 0) {
            LoadSpeedTest();
        }
        if (std::strcmp(argv[i], "--load_save") == 0 || std::strcmp(argv[i], "--all") == 0) {
            LoadSaveTest();
        }
        if (std::strcmp(argv[i], "--load_save2") == 0 || std::strcmp(argv[i], "--all") == 0) {
            LoadSaveTest2();
        }
        if (std::strcmp(argv[i], "--hops_simple") == 0 || std::strcmp(argv[i], "--all") == 0) {
            HopsTest();
        }
        if (std::strcmp(argv[i], "--hops_automorphism") == 0 || std::strcmp(argv[i], "--all") == 0) {
            HopsAutomorphismTest();
        }
        if (std::strcmp(argv[i], "--load_csv") == 0 || std::strcmp(argv[i], "--all") == 0) {
            LoadCSV();
        }
        if (std::strcmp(argv[i], "--hops_real_world") == 0 || std::strcmp(argv[i], "--all") == 0) {
            HopsRealWorldTest();
        }
        if (std::strcmp(argv[i], "--hops_simple_patterns") == 0 || std::strcmp(argv[i], "--all") == 0) {
            HopsSimplePatternsTest();
        }
        if (std::strcmp(argv[i], "--conversion") == 0 || std::strcmp(argv[i], "--all") == 0) {
            ConversionTest();
        }
        if (std::strcmp(argv[i], "--hops_time") == 0 || std::strcmp(argv[i], "--all") == 0) {
            HopsTimeTest();
        }
        if (std::strcmp(argv[i], "--hops_30s") == 0 || std::strcmp(argv[i], "--all") == 0) {
            Hops30s(threads);
        }
        if (std::strcmp(argv[i], "--hops_pattern") == 0 || std::strcmp(argv[i], "--all") == 0) {
            HopsPatternTest();
        }
        if (std::strcmp(argv[i], "--hops_parallelization") == 0 || std::strcmp(argv[i], "--all") == 0) {
            HopsParallelizationTest();
        }
        if (std::strcmp(argv[i], "--load_graphs") == 0 || std::strcmp(argv[i], "--all") == 0) {
            TestLoadGraphsFromPath("../../../cyclic_hops_cpp/data/patterns/size3", "", ".txt");
        }
        if (std::strcmp(argv[i], "--dijkstra") == 0 || std::strcmp(argv[i], "--all") == 0) {
            DijkstraTest();
        }
        if (std::strcmp(argv[i], "--load_pattern") == 0 || std::strcmp(argv[i], "--all") == 0) {
            GraphData g_patterns = GraphData<GraphStruct>("../../../cyclic_hops_cpp/data/patterns/size4/", "", "",
                                                          ".txt");
        }
        if (std::strcmp(argv[i], "--load_graph") == 0 || std::strcmp(argv[i], "--all") == 0) {
            GraphData graphs = GraphData<GraphStruct>("../../GraphData/Hops/", "", "cit-hepPh", ".bgfs");
        }
        if (std::strcmp(argv[i], "--erdos_renyi") == 0 || std::strcmp(argv[i], "--all") == 0) {
            ErdosRenyi();
        }
        if (std::strcmp(argv[i], "--dfs") == 0 || std::strcmp(argv[i], "--all") == 0) {
            DFS();
        }
        if (std::strcmp(argv[i], "--reorder") == 0 || std::strcmp(argv[i], "--all") == 0) {
            MapGraph();
        }
    }
    int x=0;
}
