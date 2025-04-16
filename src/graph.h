#ifndef GRAPH_H
#define GRAPH_H

#include <list>
#include <vector>
#include <algorithm> // Para usar std::sort
#include <set>
#include <string>

using Vertex = int;

struct Edge {
    Vertex u;
    Vertex v;
    Edge(Vertex u, Vertex v) : u(u), v(v) {}
};

class Graph {
public:
    Graph(int num_vertices);
    ~Graph();

    void add_edge(const Edge& edge);
    void print_edges() const;
    int welshPowell() const;
    int coloracaoGulosa() const;
    int dsatur() const;
    int RFL() const;
    int zykov() const;
    int lawler() const;
    int PLI() const;
    
    // Função para ler o grafo de um arquivo
    static Graph read_graph_from_file(const std::string& filename);
    
    // Função para escrever um grafo em um arquivo
    void write_graph_to_file(const std::string& filename) const;
	
private:
    int num_vertices;
    std::vector<std::vector<Vertex>> adj_list; // Lista de adjacência

    // Funções auxiliares para o Algoritmo de Zykov
    int cor(Graph& G, int q);  // Função recursiva para Zykov
    bool eh_completo();  // Verifica se o grafo é completo
    Graph contrair_vertices(int u, int v);  // Cria um grafo com contração de vértices
    Graph adicionar_aresta(int u, int v);   // Cria um grafo com adição de aresta
    
    // Função auxiliar para o algoritmo de Lawler
    bool isIndependentSet(int subset) const;
};

// Funções para gerar grafos aleatórios
double rand_double();
Graph graph_GNP(int V, double p);

// Função para executar os testes de coloração
void executarTestes(int N, double P);                // Gera grafo usando N e P
void executarTestes(const Graph& g);                 // Usa grafo já pronto


#endif // GRAPH_H
