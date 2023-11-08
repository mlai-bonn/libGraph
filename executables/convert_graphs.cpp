//
// Created by florian on 06.11.23.
//


// load all graphs from a path and convert them into bgfs format
void convert_graphs(std::string path){
    // recursively iterate over all files in the given path
    for (const auto & entry : std::filesystem::recursive_directory_iterator(path)){
        // check if the file is a graphml file
        if (entry.path().extension() == ".graphml"){
            // load the graph
            GraphStruct graph(entry.path().string());
            // save the graph in bgfs format
            graph.save(entry.path().string().substr(0, entry.path().string().size() - 8) + ".bgfs");
        }
    }
}

// main
int main(int argc, char* argv[]){
    // -p: path to the graphs
    // -j: number of threads
    // -e: file extension of the graphs
    // -f: output format
    // -o: output path
    for (int i = 0; i < argc; ++i) {
        // help
        if (strcmp(argv[i], "-h") == 0){
            std::cout << "Converts all graphs in the given path to the given output format." << std::endl;
            std::cout << "Arguments:" << std::endl;
            std::cout << "-p: path to the graphs" << std::endl;
            std::cout << "-j: number of threads (default: maximum number of threads)" << std::endl;
            std::cout << "-e: file extension of the graphs (default: .edges)" << std::endl;
            std::cout << "-f: "
            std::cout << "-o: output path (default: empty, i.e., the same as the input path)"
            return 0;
        }
    }
    return 0;
}