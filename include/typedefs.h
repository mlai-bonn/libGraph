//
// Created by Florian on 27.04.2021.
//

#ifndef HOPS_TYPEDEFS_H
#define HOPS_TYPEDEFS_H

#include <vector>
#include "Enums.h"



typedef unsigned int INDEX;
typedef INDEX Label;
typedef INDEX NodeId;
typedef std::vector<Label> Labels;
typedef std::pair<NodeId,NodeId> EDGE;
typedef std::vector<EDGE> EDGES;
typedef std::vector<NodeId> Nodes;
typedef unsigned long long int UInt64;
typedef std::vector<NodeId> PATH;
typedef std::vector<std::vector<NodeId>> PATHS;

#endif //HOPS_TYPEDEFS_H
