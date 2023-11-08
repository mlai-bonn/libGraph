//
// Created by anonymous on 12.08.21.
//

#ifndef CLOSURES_OUTERPLANARSUBGRAPH_H
#define CLOSURES_OUTERPLANARSUBGRAPH_H



class OuterplanarSubgraph {
public:
    explicit OuterplanarSubgraph(const GraphStruct & graph) : _graph(graph) {};

    //copying
    //OuterplanarSubgraph(const OuterplanarSubgraph& other){};

    virtual void generate(GraphStruct& subgraph, int seed, bool print) = 0;

    void subgraphs(int num_subgraphs, std::vector<GraphStruct>& subgraphs, int seed, const std::string& out_path = ""){
        std::mt19937_64 generator(seed);
        std::uniform_int_distribution<int> distribution(0, std::numeric_limits<int>::max());

        if (out_path.empty()) {
            subgraphs.clear();
            for (int i = 0; i < num_subgraphs; ++i) {
                subgraphs.emplace_back();
                generate(subgraphs.back(), distribution(generator), false);
            }
        }
        else{
            // delete all files in the directory that end with .bgfs
            for (const auto & entry : std::filesystem::directory_iterator(out_path)){
                if (entry.path().extension() == ".bgfs"){
                    std::filesystem::remove(entry.path());
                }
            }
            for (int i = 0; i < num_subgraphs; ++i) {
                GraphStruct subgraph;
                SaveParams saveParams = {.graphPath = out_path, .Name = "subgraph_" + std::to_string(i)};
                generate(subgraph, distribution(generator), false);
                subgraph.Save(saveParams);
            }
        }
    }

    void subgraph_extended(GraphStruct& subgraph, OuterplanarGraphData& outerplanarGraphData, int seed, bool print){
        generate(subgraph, seed, print);
        outerplanarGraphData = OuterplanarGraphData(subgraph);
    }
    void subgraphs_extended(int num_subgraphs, std::vector<OuterplanarGraphData>& subgraphs, int seed){
        std::mt19937_64 generator(seed);
        std::uniform_int_distribution<int> distribution(0, std::numeric_limits<int>::max());
        subgraphs.clear();
        for (int i = 0; i < num_subgraphs; ++i) {
            GraphStruct subgraph;
            subgraphs.emplace_back();
            subgraph_extended(subgraph, subgraphs.back(), distribution(generator), false);
        }
    }




protected:
    const GraphStruct& _graph;
};




struct OuterplanarGraphStatistics{
    OuterplanarGraphStatistics() = default;
    explicit OuterplanarGraphStatistics(const GraphStruct& outerplanarGraph);
    std::vector<std::vector<int>> ComponentSizes;
    std::vector<std::vector<int>> ComponentFaces;
    std::vector<double> ComponentNumber;
    std::vector<double> BiggestComponent;
    std::vector<double> FaceNumbers;
    void operator +=(const OuterplanarGraphStatistics& other){
        ComponentSizes.insert(ComponentSizes.end(), other.ComponentSizes.begin(),  other.ComponentSizes.end());
        ComponentNumber.insert(ComponentNumber.end(), other.ComponentNumber.begin(),  other.ComponentNumber.end());
        BiggestComponent.insert(BiggestComponent.end(), other.BiggestComponent.begin(),  other.BiggestComponent.end());
        FaceNumbers.insert(FaceNumbers.end(), other.FaceNumbers.begin(), other.FaceNumbers.end());
    };
    void getStatistics(const GraphStruct& outerplanarGraph);
    void evaluate(std::vector<std::string>& headers, std::vector<std::string>& values) const;
};

void OuterplanarGraphStatistics::getStatistics(const GraphStruct &outerplanarGraph) {
    std::vector<GraphStruct> components;
    GraphAlgorithms::GetBiconnectedComponents(outerplanarGraph, components);
    this->ComponentSizes.clear();
    this->ComponentFaces.clear();
    this->ComponentSizes.emplace_back(0);
    this->ComponentFaces.emplace_back(0);
    for (auto& component : components) {
        if (component.nodes() > 2) {
            this->ComponentSizes.back().emplace_back(component.nodes());
            this->ComponentFaces.back().emplace_back(0.0);
            GraphAlgorithms::GetBiconnectedOuterplanarFaceNum(component, this->ComponentFaces.back().back());
        }
    }
    this->ComponentNumber.emplace_back(static_cast<int>(this->ComponentSizes.back().size()));
    this->BiggestComponent.emplace_back(*std::max_element(this->ComponentSizes.back().begin(), this->ComponentSizes.back().end()));
    this->FaceNumbers.emplace_back(*std::max_element(this->ComponentFaces.back().begin(), this->ComponentFaces.back().end()));
}

OuterplanarGraphStatistics::OuterplanarGraphStatistics(const GraphStruct& outerplanarGraph) {
    getStatistics(outerplanarGraph);
}


void OuterplanarGraphStatistics::evaluate(std::vector<std::string> &headers, std::vector<std::string> &values) const {
    headers.clear();
    values.clear();
    headers = {"Avg. Component Number", "Std. Component Number", "Avg. Biggest Component", "Std. Biggest Component", "Avg. Face Number", "Std. Face Number"};
    values = {std::to_string(StaticFunctionsLib::mean(ComponentNumber)),
              std::to_string(StaticFunctionsLib::standard_deviation(ComponentNumber)),
              std::to_string(StaticFunctionsLib::mean(BiggestComponent)),
              std::to_string(StaticFunctionsLib::standard_deviation(BiggestComponent)),
              std::to_string(StaticFunctionsLib::mean(FaceNumbers)),
              std::to_string(StaticFunctionsLib::standard_deviation(FaceNumbers))};
}






#endif //CLOSURES_OUTERPLANARSUBGRAPH_H
