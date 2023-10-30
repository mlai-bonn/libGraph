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

    virtual void generate(int seed, bool print) = 0;

    GraphStruct& subgraph(int seed, bool print){
        generate(seed, print);
        return _outerPlanarSubGraph;
    }
    OuterplanarGraphData& subgraph_extended(int seed, bool print){
        generate(seed, print);
        return _outerPlanarSubGraphData;
    }


protected:
    GraphStruct _outerPlanarSubGraph;
    const GraphStruct& _graph;
    OuterplanarGraphData _outerPlanarSubGraphData;
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
    GetBiconnectedComponents(outerplanarGraph, components);
    this->ComponentSizes.clear();
    this->ComponentFaces.clear();
    this->ComponentSizes.emplace_back(0);
    this->ComponentFaces.emplace_back(0);
    for (auto& component : components) {
        if (component.nodes() > 2) {
            this->ComponentSizes.back().emplace_back(component.nodes());
            this->ComponentFaces.back().emplace_back(0.0);
            GetBiconnectedOuterplanarFaceNum(component, this->ComponentFaces.back().back());
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
