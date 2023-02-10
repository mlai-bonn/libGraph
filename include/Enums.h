//
// Created by Florian on 16.04.2021.
//

#ifndef HOPS_ENUMS_H
#define HOPS_ENUMS_H

enum ROOT_NODE_CONDITION {
    TREE_GIVEN,
    NODE_GIVEN,
    MAX_DEGREE,
    MIN_DEGREE,
    MIN_LABEL_BIG_GRAPH,
};

enum class LABEL_TYPE{
    UNLABELED,
    LABELED_SPARSE,
    LABELED_DENSE,
};

enum class RUN_TYPE{
    VECTOR,
    UNORDERED_SET,
};

enum class HOPS_TYPE{
    ESTIMATION,
    GRAPH_EDIT_DISTANCE,
};

#endif //HOPS_ENUMS_H
