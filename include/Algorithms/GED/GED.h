//
// Created by florian on 29.08.25.
//

#ifndef TUDATASETS_GED_H
#define TUDATASETS_GED_H

// include all GED related headers
#include "GEDApproximation.h"

// cpp util classes
#include "cpp_util/AscendingHeuristicEditPathComparator.h"
#include "cpp_util/AscendingOrderMyNodeComparator.h"
#include "cpp_util/CMUHouseCostFunction.h"
#include "cpp_util/Constants.h"
#include "cpp_util/DescendingHeuristicEditPathComparator.h"
#include "cpp_util/Edge.h"
#include "cpp_util/EditPath.h"
#include "cpp_util/Graph.h"
#include "cpp_util/GraphCollection.h"
#include "cpp_util/GraphComponent.h"
#include "cpp_util/GRECCostFunction.h"
#include "cpp_util/ICostFunction.h"
#include "cpp_util/IEdgeHandler.h"
#include "cpp_util/INodeHandler.h"
#include "cpp_util/Matrix.h"
#include "cpp_util/MatrixElement.h"
#include "cpp_util/MutagenCostFunction.h"
#include "cpp_util/MyTree.h"
#include "cpp_util/Node.h"
#include "cpp_util/UniversalEdgeHandler.h"
#include "cpp_util/UnlabeledCostFunction.h"

// cpp
#include "cpp_algorithms/AscendingOrderNonStaticNodeComparator.h"
#include "cpp_algorithms/NonStaticTree.h"
#include "cpp_algorithms/Constants.h"
#include "cpp_algorithms/GEDDFS.h"
#include "cpp_algorithms/GEDDFSThread.h"
#include "cpp_algorithms/GEDMultiThread.h"
#include "cpp_algorithms/GlobalVar.h"
#include "cpp_algorithms/HeuristicGEDPartialMatcher.h"
#include "cpp_algorithms/MatrixGenerator.h"
#include "cpp_algorithms/Munkres.h"
#include "cpp_algorithms/MunkresRec.h"
#include "cpp_algorithms/NonStaticNode.h"
#include "cpp_algorithms/NonStaticTree.h"

#endif //TUDATASETS_GED_H