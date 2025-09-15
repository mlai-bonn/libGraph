//
// Created by Florian on 22.04.2021.
//

#ifndef HOPS_STATICFUNCTIONS_H
#define HOPS_STATICFUNCTIONS_H

class GraphStaticFunctions {
public:
    static void ConvertWikiTopcats(const std::string& graphPath, const std::string& categoryPath){
        GraphStruct graph = GraphStruct(graphPath, true, true, categoryPath);
        graph.Save({std::filesystem::path(graphPath).parent_path().string() + "/", "", GraphFormat::BINARY, true});
    }
};


#endif //HOPS_STATICFUNCTIONS_H
