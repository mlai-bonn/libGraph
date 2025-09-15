//
// Created by Florian on 15.04.2021.
//

#ifndef HOPS_LOAD_H
#define HOPS_LOAD_H

#include <string>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstring>
//    labels_from_gml("D:\\EigeneDokumente\\Forschung\\Code\\cyclic_hops\\data\\graphs\\soc-pokec\\soc-pokec.gml", "D:\\EigeneDokumente\\Forschung\\Code\\cyclic_hops\\data\\graphs\\soc-pokec\\soc-pokec.labels");

static void removeCharsFromString(std::string &str,const char* charsToRemove ) {
    for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) {
        str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
    }
}

static void labels_from_gml(std::string input_path, std::string output_path){
    std::ifstream infile(input_path);
    std::string line;
    size_t counter = 0;
    std::ofstream out(output_path);
    std::string key;
    std::string value;
    size_t label;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        iss >> key >> value;
        if (key == "label"){
            removeCharsFromString(value, "\"");
            std::stringstream sstream(value);
            sstream >> label;
            out << counter << " " << label << std::endl;
            ++counter;
        }
    }
    out.close();
}




#endif //HOPS_LOAD_H
