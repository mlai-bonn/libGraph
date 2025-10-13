//
// Created by Florian on 16.04.2021.
//

#ifndef LIBGRAPH_ENUMS_H
#define LIBGRAPH_ENUMS_H

#include <iostream>

enum ROOT_NODE_CONDITION {
    TREE_GIVEN,
    NODE_GIVEN,
    MAX_DEGREE,
    MIN_DEGREE,
    MIN_LABEL_BIG_GRAPH,
};

inline std::ostream& operator<<(std::ostream& os, const ROOT_NODE_CONDITION& condition) {
    switch (condition) {
        case ROOT_NODE_CONDITION::TREE_GIVEN:
            return os << "TREE_GIVEN";
        case ROOT_NODE_CONDITION::NODE_GIVEN:
            return os << "NODE_GIVEN";
        case ROOT_NODE_CONDITION::MAX_DEGREE:
            return os << "MAX_DEGREE";
        case ROOT_NODE_CONDITION::MIN_DEGREE:
            return os << "MIN_DEGREE";
        case ROOT_NODE_CONDITION::MIN_LABEL_BIG_GRAPH:
            return os << "MIN_LABEL_BIG_GRAPH";
    }
    return os;
}

enum class LABEL_TYPE{
    UNLABELED,
    LABELED_SPARSE,
    LABELED_DENSE,
};

inline std::ostream& operator<<(std::ostream& os, const LABEL_TYPE& type) {
    switch (type) {
        case LABEL_TYPE::UNLABELED:
            return os << "UNLABELED";
        case LABEL_TYPE::LABELED_SPARSE:
            return os << "LABELED_SPARSE";
        case LABEL_TYPE::LABELED_DENSE:
            return os << "LABELED_DENSE";
    }
    return os;
}

enum class RUN_TYPE{
    VECTOR,
    UNORDERED_SET,
};

inline std::ostream& operator<<(std::ostream& os, const RUN_TYPE& type) {
    switch (type) {
        case RUN_TYPE::VECTOR:
            return os << "VECTOR";
        case RUN_TYPE::UNORDERED_SET:
            return os << "UNORDERED_SET";
    }
    return os;
}

enum class HOPS_TYPE{
    ESTIMATION,
    GRAPH_EDIT_DISTANCE,
};

inline std::ostream& operator<<(std::ostream& os, const HOPS_TYPE& type) {
    switch (type) {
        case HOPS_TYPE::ESTIMATION:
            return os << "ESTIMATION";
        case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
            return os << "GRAPH_EDIT_DISTANCE";
    }
    return os;
}

#endif //LIBGRAPH_ENUMS_H
