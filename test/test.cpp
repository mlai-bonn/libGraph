
#include "load_save.h"
#include "graph_functions.h"
#include "hops.h"

int main() {
    //GraphData g_patterns = GraphData<GraphStruct>("../../../GraphData/Hops/patterns/size3.bgf");
    //LoadSpeedTest();
    //LoadSaveTest();
    //LoadSaveTest2();
    //HopsTest();
    //HopsAutomorphismTest();
    //LoadCSV();
    //HopsRealWorldTest();
    //HopsSimplePatternsTest();
    //ConversionTest();
    HopsTimeTest();
    //HopsPatternTest();
    //HopsParallelizationTest();
    //TestLoadGraphsFromPath("../../../cyclic_hops_cpp/data/patterns/size3","", ".txt");
    //DijkstraTest();
    //GraphData g_patterns = GraphData("../../../cyclic_hops_cpp/data/patterns/size4/", "", "", ".txt");
    GraphData graphs = GraphData<GraphStruct>("../../GraphData/Hops/", "", "cit-hepPh", ".bgfs");
    int x=0;
}
