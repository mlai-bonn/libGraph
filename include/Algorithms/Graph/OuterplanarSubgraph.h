//
// Created by anonymous on 12.08.21.
//

#ifndef CLOSURES_OUTERPLANARSUBGRAPH_H
#define CLOSURES_OUTERPLANARSUBGRAPH_H

class OuterplanarSubgraph {
public:
    explicit OuterplanarSubgraph(const GraphStruct & graph);

    //copying
    OuterplanarSubgraph(const OuterplanarSubgraph& other);

    virtual void generate(GraphStruct& subgraph, std::mt19937_64& gen, bool print) = 0;

protected:
    std::mt19937_64 _gen;
    const GraphStruct _graph;
    bool print = false;
};

OuterplanarSubgraph::OuterplanarSubgraph(const GraphStruct& graph) : _graph(graph){}

OuterplanarSubgraph::OuterplanarSubgraph(const OuterplanarSubgraph &other) : _graph(other._graph), _gen(other._gen), print(other.print) {
}

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
    std::vector<std::vector<NodeId>> components;
    GetBiconnectedComponents(outerplanarGraph, components);
    this->ComponentSizes.clear();
    this->ComponentFaces.clear();
    this->ComponentSizes.emplace_back(0);
    this->ComponentFaces.emplace_back(0);
    for (const auto& component : components) {
        if (component.size() > 2) {
            this->ComponentSizes.back().emplace_back(component.size());
            this->ComponentFaces.back().emplace_back(0.0);
            GraphFunctions::GetBiconnectedOuterplanarFaceNum(component, this->ComponentFaces.back().back());
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
