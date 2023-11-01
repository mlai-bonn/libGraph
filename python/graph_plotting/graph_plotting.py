import matplotlib.pyplot as plt
import networkx as nx


def main():
    # load graph from path
    path = "../../Google_tests/tests/out/"
    # get erdos renyi graph
    er_graph = nx.read_edgelist(path + "erdos_renyi.edges")
    er_outerplanar = nx.read_edgelist(path + "erdos_renyi_outerplanar_subgraph.edges")

    pos = {'0' : (0,0),
           '1' : (2,0),
           '2' : (1,0),
           '3' : (0,1),
           '4' : (2,1),
           '5' : (1,1),
           '6' : (0,2),
           '7' : (2,2),
           '8' : (1,2),
           '9' : (1,0)}
    pos = nx.spring_layout(er_graph, seed=203745)

    nx.draw_networkx_nodes(er_graph, pos, node_size=10)
    nx.draw_networkx_edges(er_graph, pos, width=0.1)
    plt.show()
    nx.draw_networkx_nodes(er_outerplanar, pos, node_size=10)
    nx.draw_networkx_edges(er_outerplanar, pos, width=0.1)
    plt.show()

# run main
if __name__ == "__main__":
    main()
