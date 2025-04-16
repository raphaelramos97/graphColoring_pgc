#include <iostream>
#include <iomanip>
#include "graph.h"

int main() {
    srand(time(0));

    // Grupo de testes 1: grafo aleatório gerado pro graph_GNP()
    int valoresN[] = {5,10,15};
    double valoresP[] = {0.2,0.4,0.6,0.8};

    for (int N : valoresN) {
        for (double P : valoresP) {
			Graph g = graph_GNP(N, P); 
			//g.write_graph_to_file("output_graph.col");
            executarTestes(g);
        }
    }

    // Grupo de testes 2: instâncias DIMACS 
    Graph graph = Graph::read_graph_from_file("myciel3.col");
    executarTestes(graph); 
    

    return 0;
}
