//
// Created by Florian on 12.04.2021.
//

#ifndef CLOSURES_BASEOPERATOR_H
#define CLOSURES_BASEOPERATOR_H

#include <set>
#include <utility>
#include "../io/StaticFunctions.h"

struct ClosureParameters {
public:
    std::set<NodeId> input_set  = std::set<NodeId>();
    int threshold = std::numeric_limits<int>::max();
    bool output_is_all = false;
    bool output_intersects_forbidden = false;
    NodeId element_to_add = -1;
    bool onlyPreClosure = false;
    int incrementalCount = -1;
    double timeConstraint = -1;
    double lower_bound = 0.75;
    int iteration_number = 0;
    int target_set_size = 0;
    int preclosure_depth = 0;
    bool detailed_analysis = false;
    std::map<int, std::vector<NodeId>> preClosureDepthToElements;

    std::set<NodeId> closed_set = std::set<NodeId>();
    std::set<NodeId> added_elements = std::set<NodeId>();
    const std::set<NodeId>* forbidden_elements = nullptr;





    void print()
    {
        std::cout << "Size: " << closed_set.size() << std::endl;
        std::cout << "Input Set: " << StaticFunctionsLib::print<std::set<NodeId>, NodeId>(this->input_set) << std::endl;
        std::cout << "Closed Set: " << StaticFunctionsLib::print<std::set<NodeId>, NodeId>(this->closed_set) << std::endl;
        std::cout << "Added Elements: " << StaticFunctionsLib::print<std::set<NodeId>, NodeId>(this->added_elements) << std::endl;
    };
    void reset()
    {
        this->output_intersects_forbidden = false;
        this->output_is_all = false;
        lower_bound = 0.75;
        target_set_size = 0;
        preclosure_depth = 0;
    };


    void clear(){
        this->output_intersects_forbidden = false;
        this->output_is_all = false;
        this->onlyPreClosure = false;
        this->incrementalCount = -1;
        this->timeConstraint = -1;
        this->closed_set.clear();
        this->added_elements.clear();
        this->preClosureDepthToElements.clear();
        lower_bound = 0.75;
        target_set_size = 0;
        iteration_number = 0;
    };
};

template <typename T>
class BaseOperator {
public:
    explicit BaseOperator(std::string name) : _name(std::move(name)){}

    BaseOperator();;
    std::string getName(){return _name;};

    virtual void closure(ClosureParameters& ClosureOutput) = 0;


protected:
    [[nodiscard]] bool IsForbidden(NodeId Id, const std::vector<std::set<NodeId>*>& forbidden_elements) const;

private:
    std::string _name;

};

template<typename T>
bool BaseOperator<T>::IsForbidden(NodeId Id, const std::vector<std::set<NodeId>*>& forbidden_elements) const {
    for (std::set<NodeId>* nodeSet : forbidden_elements) {
        if (nodeSet->find(Id) != nodeSet->end()){
            return true;
        }
    }
    return false;
}

template<typename T>
BaseOperator<T>::BaseOperator() = default;


#endif //CLOSURES_BASEOPERATOR_H
