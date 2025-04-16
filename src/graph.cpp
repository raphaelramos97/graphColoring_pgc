#include "graph.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>
#include <queue>
#include <set>
#include <climits> 
#include <cassert>
#include <cstdlib>  // Para rand() e RAND_MAX
#include <ctime>    // Para srand()
#include <cmath>    // Para rand_double()
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <glpk.h>


// Construtor do grafo
Graph::Graph(int num_vertices) : num_vertices(num_vertices), adj_list(num_vertices) {}

// Destrutor do grafo
Graph::~Graph() {}

// Adiciona uma aresta ao grafo
void Graph::add_edge(const Edge& edge) {
    if (edge.u >= 0 && edge.u < num_vertices && edge.v >= 0 && edge.v < num_vertices) {
        adj_list[edge.u].push_back(edge.v);
        adj_list[edge.v].push_back(edge.u); // Se o grafo for não dirigido
    } else {
        std::cerr << "Erro: Vértice fora dos limites!" << std::endl;
    }
}


// Imprime todas as arestas do grafo
void Graph::print_edges() const {
    for (Vertex u = 0; u < num_vertices; ++u) {
        for (Vertex v : adj_list[u]) {
            if (u < v) { // Para evitar a impressão de arestas duplicadas em grafos não dirigidos
                std::cout << "(" << u << ", " << v << ")" << std::endl;
            }
        }
    }
}

// Verifica se o grafo é completo
bool Graph::eh_completo() {
    for (int i = 0; i < num_vertices; ++i) {
        if (adj_list[i].size() != num_vertices - 1) {
            return false;
        }
    }
    return true;
}

// Contrai dois vértices u e v em um grafo
Graph Graph::contrair_vertices(int u, int v) {
    Graph G_contraido(num_vertices - 1);  // O novo grafo terá um vértice a menos
    std::vector<int> novo_mapeamento(num_vertices, -1);  // Mapeia os vértices antigos para os novos
    std::set<std::pair<int, int>> arestas_adicionadas;  // Para evitar duplicação de arestas

    // Cria um novo mapeamento de vértices, onde u e v são mesclados
    int novo_vertice = 0;
    for (int i = 0; i < num_vertices; ++i) {
        if (i != u && i != v) {
            novo_mapeamento[i] = novo_vertice++;
        }
    }

    // Insere as arestas no novo grafo, com os vértices remapeados
    for (int i = 0; i < num_vertices; ++i) {
        if (i == u || i == v) {
            continue;  // O vértice u ou v será substituído pelo vértice mesclado
        }
        
        for (int j : adj_list[i]) {
            int v1 = novo_mapeamento[i];
            int v2 = (j == u || j == v) ? novo_vertice : novo_mapeamento[j];
            
            // Armazena a aresta de forma ordenada para evitar duplicação
            if (v1 != v2) {
                auto aresta = std::minmax(v1, v2);  // Garante que as arestas sejam ordenadas
                if (arestas_adicionadas.insert(aresta).second) {
                    G_contraido.add_edge(Edge(aresta.first, aresta.second));
                }
            }
        }
    }

    return G_contraido;
}


// Adiciona uam nova aresta (u, v) no grafo
Graph Graph::adicionar_aresta(int u, int v) {
    Graph G_adicionado = *this;
    G_adicionado.add_edge(Edge(u, v));
    return G_adicionado;
}

// Gera um número aleatório entre 0 e 1
double rand_double() {
    return static_cast<double>(rand()) / RAND_MAX;
}

// Gerar o grafo GNP
Graph graph_GNP(int V, double p) {
    assert(V > 0);
    assert(p >= 0 && p <= 1);

    Graph G(V);  // Inicializa o grafo com V vértices

    // Aresta uv existe com probabilidade p
    for (int u = 0; u < V; u++) {
        for (int v = u + 1; v < V; v++) {
            if (rand_double() <= p) {
                G.add_edge(Edge(u, v));
            }
        }
    }

    return G;
}

// Função para verificação de conjuntos independentes
bool Graph::isIndependentSet(int subset) const {
    for (int u = 0; u < num_vertices; ++u) {
        if (subset & (1 << u)) { // Verifica se u está no subconjunto
            for (int v : adj_list[u]) {
                if (subset & (1 << v)) { // Se v também estiver no subconjunto, não é independente
                    return false;
                }
            }
        }
    }
    return true;
}

// Função para ler o grafo de um arquivo
Graph Graph::read_graph_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erro: Não foi possível abrir o arquivo!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    int num_vertices = 0, num_edges = 0;

    // Primeira leitura para obter o número de vértices
    while (std::getline(file, line)) {
        // Ignora comentários
        if (line[0] == 'c') {
            continue;
        }

        // Linha com o número de vértices e arestas
        if (line[0] == 'p') {
            std::istringstream iss(line);
            std::string ignore, type;
            iss >> ignore >> type >> num_vertices >> num_edges;
            break; // Após encontrar essa linha, podemos parar de procurar
        }
    }

    // Cria o grafo com o número de vértices lido
    Graph graph(num_vertices);

    // Segunda leitura para obter as arestas
    while (std::getline(file, line)) {
        if (line[0] == 'e') {
            std::istringstream iss(line);
            char e;
            int u, v;
            iss >> e >> u >> v;
            u--; v--; // Ajuste para indexação começando em 0

            Edge edge(u, v);
            graph.add_edge(edge); // Usa a função de adicionar aresta no grafo
        }
    }

    file.close();
    return graph;
}

// Função para escrever um grafo de um arquivo
void Graph::write_graph_to_file(const std::string& filename) const {
    std::ofstream file(filename);  // Abre o arquivo para escrita

    if (!file.is_open()) {
        std::cerr << "Erro: Não foi possível abrir o arquivo para escrita!" << std::endl;
        return;
    }

    // Escreve comentários no início do arquivo
    file << "c Arquivo gerado automaticamente\n";
    file << "c Este grafo tem " << num_vertices << " vértices\n";

    // Conta o número de arestas
    int num_edges = 0;
    for (int u = 0; u < num_vertices; ++u) {
        num_edges += adj_list[u].size();
    }
    num_edges /= 2;  // Cada aresta é contada duas vezes (uma para cada vértice)

    // Escreve a linha 'p edge' com o número de vértices e arestas
    file << "p EDGE " << num_vertices << " " << num_edges << "\n";

    // Escreve as arestas no formato 'e u v'
    for (int u = 0; u < num_vertices; ++u) {
        for (int v : adj_list[u]) {
            if (u < v) {  // Para evitar duplicar as arestas
                file << "e " << (u + 1) << " " << (v + 1) << "\n";  // Aumenta os índices em 1 para seguir o formato do arquivo (vértices começam em 1)
            }
        }
    }

    file.close();  // Fecha o arquivo
}



// Função para executar os testes em um grafo gerado pela função graph_GNP
void executarTestes(int N, double P) {
    Graph g = graph_GNP(N, P); 
    executarTestes(g);
}

// Função para executar os testes em um grafo pronto
void executarTestes(const Graph& g) {

		
    // Medindo tempo e número de cores do Algoritmo Guloso
    auto start = std::chrono::high_resolution_clock::now();
    int coresGuloso = g.coloracaoGulosa();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracaoGuloso = (end - start);

    // Medindo tempo e número de cores do Algoritmo Welsh-Powell
    start = std::chrono::high_resolution_clock::now();
    int coresWelshPowell = g.welshPowell();
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracaoWelshPowell = (end - start);

    // Medindo tempo e número de cores do Algoritmo DSATUR
    start = std::chrono::high_resolution_clock::now();
    int coresDSATUR = g.dsatur();
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracaoDSATUR = (end - start);

    // Medindo tempo e número de cores do Algoritmo RFL
    start = std::chrono::high_resolution_clock::now();
    int coresRFL = g.RFL();
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracaoRFL = (end - start);

    // Medindo tempo e número de cores do Algoritmo Zykov
    start = std::chrono::high_resolution_clock::now();
    int coresZykov = g.zykov();
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracaoZykov = (end - start);

    // Medindo tempo e número de cores do Algoritmo de Lawler
	start = std::chrono::high_resolution_clock::now();
	int coresLawler = g.lawler();
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duracaoLawler = (end - start);
	
	// Medindo tempo e número de cores do Algoritmo PLI
	start = std::chrono::high_resolution_clock::now();
	int coresPLI = g.PLI();
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duracaoPLI = (end - start);

    // Exibindo os resultados
    //std::cout << "Resultados para N = " << N << ", P = " << P << ":\n";
    std::cout << "------------------------------------------------------------------\n";
    std::cout << std::setw(20) << "Algoritmo" 
              << std::setw(20) << "Qtde. Cores"
              << std::setw(20) << "Tempo (s)" << std::endl;
    
    std::cout << std::setw(20) << "Guloso" 
              << std::setw(20) << coresGuloso
              << std::setw(20) << duracaoGuloso.count() << std::endl;

    std::cout << std::setw(20) << "Welsh-Powell" 
              << std::setw(20) << coresWelshPowell
              << std::setw(20) << duracaoWelshPowell.count() << std::endl;

    std::cout << std::setw(20) << "DSATUR" 
              << std::setw(20) << coresDSATUR
              << std::setw(20) << duracaoDSATUR.count() << std::endl;

    std::cout << std::setw(20) << "RFL" 
              << std::setw(20) << coresRFL
              << std::setw(20) << duracaoRFL.count() << std::endl;
              
    std::cout << std::setw(20) << "Zykov" 
              << std::setw(20) << coresZykov
              << std::setw(20) << duracaoZykov.count() << std::endl;
              
    std::cout << std::setw(20) << "Lawler" 
			  << std::setw(20) << coresLawler
              << std::setw(20) << duracaoLawler.count() << std::endl;
              
    std::cout << std::setw(20) << "PLI" 
			  << std::setw(20) << coresPLI
              << std::setw(20) << duracaoPLI.count() << std::endl;

    std::cout << "------------------------------------------------------------------\n";
}

// ---------------------------- Coloração Gulosa Sequencial ----------------------------

int Graph::coloracaoGulosa() const {

    // vetor C: a posição C[i] indica a cor do vértice i
    std::vector<int> C(num_vertices, -1);
    C[0] = 0;

    // vetor disp: indica se uma cor k está disponível para uso (disp[k] == true)
    std::vector<bool> disp(num_vertices, true);

    int max_cor = 0; // variável para acompanhar a maior cor usada

    // percorremos sequencialmente os vértices de G
    for (int u = 1; u < num_vertices; ++u) {

        // para cada vértice u de V(G), percorremos N(u) e marcamos como false as cores já utilizadas
        for (Vertex vizinho : adj_list[u]) {
            if (C[vizinho] != -1) {
                disp[C[vizinho]] = false;
            }
        }

        // capturamos a menor cor disponível para uso (variável cor)
        int cor;
        for (cor = 0; cor < num_vertices; ++cor) {
            if (disp[cor]) {
                break;
            }
        }

        // colorimos o vértice u com a menor cor disponível
        C[u] = cor;

        // Atualiza a cor máxima usada
        if (cor > max_cor) {
            max_cor = cor;
        }

        // repassamos por N(u) para desmarcar os valores alterados no vetor disp (reset)
        for (Vertex vizinho : adj_list[u]) {
            if (C[vizinho] != -1) {
                disp[C[vizinho]] = true;
            }
        }
    }


    // O número cromático é a maior cor usada + 1
    return max_cor + 1;
}


// ---------------------------- Algoritmo de Welsh e Powell ----------------------------

int Graph::welshPowell() const {

    // Vetor para armazenar os graus dos vértices
    std::vector<std::pair<Vertex, int>> grau(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
        grau[i] = {i, adj_list[i].size()};
    }

    // Ordena os vértices por grau decrescente
    std::sort(grau.begin(), grau.end(), [](const std::pair<Vertex, int>& a, const std::pair<Vertex, int>& b) {
        return b.second < a.second;
    });

    // Vetor para armazenar as cores dos vértices
    std::vector<int> C(num_vertices, -1);

    // Cor inicial
    int k = 0;

    // Percorre V(G)
    for (const auto& [u, _] : grau) {

        // Pega um vértice u ainda não colorido e atribui a cor k
        if (C[u] == -1) {
            C[u] = k;

            // Pega um outro vértice v ainda não colorido
            for (const auto& [v, _] : grau) {
                if (C[v] == -1) {

                    // Verifica se algum vizinho de v já possui a cor k
                    bool disp = true;
                    for (Vertex vizinho : adj_list[v]) {
                        if (C[vizinho] == k) {

                            // Se possuir, vamos para a próxima cor
                            disp = false;
                            break;
                        }
                    }

                    // Se não possui, atribuimos a cor k ao vértice V
                    if (disp) {
                        C[v] = k;
                    }
                }
            }

            // próxima cor
            k++;
        }
    }

    // O número cromático é o valor de k (última cor utilizada)
    return k;
}



//---------------------------- Algoritmo de coloração DSATUR ----------------------------

int Graph::dsatur() const {

    // Vetor para armazenar as cores dos vértices
    std::vector<int> C(num_vertices, -1);

    // Vetor para armazenar o grau de saturação dos vértices
    std::vector<int> saturacao(num_vertices, 0);

    // Vetor para armazenar o grau dos vértices
    std::vector<int> grau(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
        grau[i] = adj_list[i].size();
    }

    // Iniciamos atribuindo a primeira cor ao vértice de maior grau
    int maiorGrau = std::max_element(grau.begin(), grau.end()) - grau.begin();
    C[maiorGrau] = 0;

    // Variável para rastrear a maior cor utilizada
    int max_cor = 0;

    // Atualizamos o grau de saturação dos vizinhos desse primeiro vértice colorido
    for (Vertex vizinho : adj_list[maiorGrau]) {
        saturacao[vizinho]++;
    }

    for (int i = 1; i < num_vertices; ++i) {

        // inicia em -1 o vértice com maior saturação e o seu grau de saturação
        int maiorSaturVertice = -1;
        int maiorSatur = -1;

        // Para cada vértice ainda não colorido, vamos procurar o que tem o maior grau de saturação
        for (int j = 0; j < num_vertices; ++j) {
            if (C[j] == -1) {
                // Se o vértice j tiver o maior grau de saturação encontrado até o momento, ele é escolhido
                if (saturacao[j] > maiorSatur) {
                    maiorSatur = saturacao[j];
                    maiorSaturVertice = j;

                // Se o grau de saturação de j for igual ao maior já encontrado, desempatamos pela escolha do grau
                } else if (saturacao[j] == maiorSatur) {
                    if (grau[j] > grau[maiorSaturVertice]) {
                        maiorSaturVertice = j;
                    }
                }
            }
        }

        // Marcamos todos os vizinhos do vértice com o maior grau de saturação como indisponíveis
        std::vector<bool> disp(num_vertices, true);
        for (Vertex vizinho : adj_list[maiorSaturVertice]) {
            if (C[vizinho] != -1) {
                disp[C[vizinho]] = false;
            }
        }

        // Capturamos a menor cor disponível para uso (variável cor)
        int cor;
        for (cor = 0; cor < num_vertices; ++cor) {
            if (disp[cor]) {
                break;
            }
        }

        // Atribuímos a menor cor disponível ao vértice com maior grau de saturação
        C[maiorSaturVertice] = cor;

        // Atualizamos a maior cor usada
        if (cor > max_cor) {
            max_cor = cor;
        }

        // Atualizamos o grau de saturação dos vizinhos (ainda não coloridos) desse vértice recém-colorido
        for (Vertex vizinho : adj_list[maiorSaturVertice]) {
            if (C[vizinho] == -1) {
                saturacao[vizinho]++;
            }
        }
    }

    // O número cromático é a maior cor usada + 1
    return max_cor + 1;
}

// ---------------------------------- Algoritmo RFL  ----------------------------------

int Graph::RFL() const {
    std::vector<std::set<int>> C;  // Partição final em conjuntos independentes
    std::set<int> V_prime;         // Conjunto de vértices não processados
    std::set<int> A;               // Conjunto de vértices adjacentes a vértices já processados

    // Inicializa V' como o conjunto de todos os vértices
    for (int i = 0; i < num_vertices; ++i) {
        V_prime.insert(i);
    }

    int k = 0; // Contador de subconjuntos independentes

    // Enquanto ainda houver vértices não processados
    while (!V_prime.empty()) {
        std::set<int> C_k;  // Conjunto independente atual
        std::set<int> A;    // Conjunto de vizinhos dos vértices processados nesta rodada

        // Processa o subconjunto independente até que não seja possível adicionar mais vértices a C_k
        std::set<int> to_remove;  // Conjunto para armazenar vértices que serão removidos de V'

        for (int v : V_prime) {
            // Verifica se algum vizinho de v já está no conjunto independente atual
            bool pode_adicionar = true;
            for (int vizinho : adj_list[v]) {
                if (C_k.find(vizinho) != C_k.end()) {
                    pode_adicionar = false;
                    break;
                }
            }

            // Se o vértice v não possui vizinhos em C_k, pode ser adicionado
            if (pode_adicionar) {
                C_k.insert(v);
                // Adiciona os vizinhos de v ao conjunto A
                for (int vizinho : adj_list[v]) {
                    A.insert(vizinho);
                }
                to_remove.insert(v);  // Marca v para remoção de V_prime
            }
        }

        // Remove os vértices processados nesta rodada de V_prime
        for (int v : to_remove) {
            V_prime.erase(v);
        }

        // Adiciona o conjunto C_k à partição C
        C.push_back(C_k);
        k++; // Incrementa o contador de conjuntos independentes
    }

    // O número cromático é o número de conjuntos independentes (partições)
    return k;
}

// ---------------------------- Algoritmo de Zykov ----------------------------

// Função recursiva do algoritmo de Zykov
int Graph::cor(Graph& G, int q){
    int n = G.num_vertices;

    // Se o grafo é completo, o número cromático é o número de vértices
    if (G.eh_completo()) {
        return n;
    } else {
        // Encontrar dois vértices não adjacentes
        int u = -1, v = -1;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (std::find(G.adj_list[i].begin(), G.adj_list[i].end(), j) == G.adj_list[i].end()) {
                    u = i;
                    v = j;
                    break;
                }
            }
            if (u != -1 && v != -1) break;
        }

        // Realizar a recursão para o grafo contraído e o grafo com a aresta adicionada
        Graph G_contraido = G.contrair_vertices(u, v);
        int q1 = cor(G_contraido, q);

        Graph G_adicionado = G.adicionar_aresta(u, v);
        int q2 = cor(G_adicionado, q);

        // Retornar o mínimo dos dois resultados
        return std::min(q1, q2);
    }
}

// Função principal do algoritmo de Zykov
int Graph::zykov() const {
    Graph copia = *this;
    return copia.cor(copia, num_vertices);
}

/*
 * int Graph::zykov() {
 *    return cor(*this, num_vertices);
 * } */


// ---------------------------- Algoritmo de Lawler ----------------------------

int Graph::lawler() const {
    int n = num_vertices;
    int tamanho_X = 1 << n; // 2^n
    std::vector<int> X(tamanho_X, INT_MAX);
    X[0] = 0; // Inicializa o primeiro elemento com 0

    // Percorre todos os subconjuntos possíveis
    for (int i = 1; i < tamanho_X; ++i) {
        // Para cada subconjunto i, encontra os subconjuntos independentes maximais
        for (int j = i; j > 0; j = (j - 1) & i) {
            if (isIndependentSet(j)) {
                // Atualiza o valor de X[i] usando o conjunto independente
                X[i] = std::min(X[i], X[i ^ j] + 1);
            }
        }
    }

    // O resultado final será o valor de X que corresponde ao conjunto completo de vértices
    return X[tamanho_X - 1];
}

// ---------------------------- Algoritmo PLI ----------------------------

int Graph::PLI() const {
    glp_prob *lp = glp_create_prob();
    glp_set_prob_name(lp, "Graph Coloring");
    glp_set_obj_dir(lp, GLP_MIN); // Minimizar número de cores

    int numY = num_vertices * num_vertices;  // Variáveis y_vc
    int numX = num_vertices;                 // Variáveis x_c
    int numVars = numY + numX;

    glp_add_cols(lp, numVars);

    // Definir variáveis y_vc (binárias)
    for (int v = 0; v < num_vertices; ++v) {
        for (int c = 0; c < num_vertices; ++c) {
            int varIndex = v * num_vertices + c + 1;
            glp_set_col_kind(lp, varIndex, GLP_BV);
        }
    }

    // Definir variáveis x_c (binárias)
    for (int c = 0; c < num_vertices; ++c) {
        int varIndex = numY + c + 1;
        glp_set_col_kind(lp, varIndex, GLP_BV);
    }

    // Restrição (2): Cada vértice recebe exatamente uma cor
    glp_add_rows(lp, num_vertices);
    for (int v = 0; v < num_vertices; ++v) {
        int rowIndex = v + 1;
        glp_set_row_bnds(lp, rowIndex, GLP_FX, 1.0, 1.0);

        int ind[num_vertices + 1];
        double val[num_vertices + 1];

        for (int c = 0; c < num_vertices; ++c) {
            ind[c + 1] = v * num_vertices + c + 1;
            val[c + 1] = 1.0;
        }
        glp_set_mat_row(lp, rowIndex, num_vertices, ind, val);
    }

    // Restrição (3): Vértices adjacentes não podem ter a mesma cor
    int rowCount = num_vertices;
    for (int u = 0; u < num_vertices; ++u) {
        for (int v : adj_list[u]) {
            if (u < v) { // Evita duplicação
                for (int c = 0; c < num_vertices; ++c) {
                    rowCount++;
                    glp_add_rows(lp, 1);
                    glp_set_row_bnds(lp, rowCount, GLP_UP, 0.0, 1.0);

                    int ind[3] = {0, u * num_vertices + c + 1, v * num_vertices + c + 1};
                    double val[3] = {0, 1.0, 1.0};

                    glp_set_mat_row(lp, rowCount, 2, ind, val);
                }
            }
        }
    }

    // Restrição (4): Se um vértice v recebe a cor c, então x_c deve ser 1
    for (int v = 0; v < num_vertices; ++v) {
        for (int c = 0; c < num_vertices; ++c) {
            rowCount++;
            glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, rowCount, GLP_UP, 0.0, 0.0);

            int ind[3] = {0, v * num_vertices + c + 1, numY + c + 1};
            double val[3] = {0, 1.0, -1.0};

            glp_set_mat_row(lp, rowCount, 2, ind, val);
        }
    }

    // Função objetivo: minimizar o número de cores usadas (soma de x_c)
    glp_set_obj_coef(lp, 0, 0.0);
    for (int c = 0; c < num_vertices; ++c) {
        glp_set_obj_coef(lp, numY + c + 1, 1.0);
    }


    // Forçar o uso do Primal Simplex antes da otimização inteira
    glp_smcp simplex_params;
    glp_init_smcp(&simplex_params);
    simplex_params.msg_lev = GLP_MSG_OFF;
    simplex_params.meth = GLP_PRIMAL; // Força o método primal simplex
    glp_adv_basis(lp, 0); // Ajusta a base inicial
    glp_simplex(lp, &simplex_params);

    // Resolver diretamente como MIP
    glp_iocp intopt_params;
    glp_init_iocp(&intopt_params);
    intopt_params.msg_lev = GLP_MSG_OFF;
    intopt_params.presolve = GLP_ON; // Ativar pré-solução
    intopt_params.binarize = GLP_ON; // Melhorar tratamento das variáveis binárias
    glp_intopt(lp, &intopt_params); // Resolve o PLI

    // Extração da solução ótima
    int chromaticNumber = 0;
    for (int c = 0; c < num_vertices; ++c) {
        if (glp_mip_col_val(lp, numY + c + 1) > 0.5) {
            chromaticNumber++;
        }
    }

    glp_delete_prob(lp);
    return chromaticNumber;
}
