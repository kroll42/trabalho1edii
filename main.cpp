#include <stdio.h>
#include <stdlib.h>
#include "grafo.h"

int main() {
    int representation;
    int startVertex;

    printf("Escolha a representacao do grafo:\n");
    printf("0 - Lista de Adjacencia\n");
    printf("1 - Matriz de Adjacencia\n");
    scanf("%d", &representation);

    printf("Digite o vertice inicial para DFS e BFS: ");
    scanf("%d", &startVertex);

    const char* statsFilename = "estatisticas.txt";
    const char* bfsFilename = "resultado_BFS.txt";
    const char* dfsFilename = "resultado_DFS.txt";
    const char* inputFilename = "grafo_2.txt"; // Substitua pelo caminho do seu arquivo de entrada

    Graph* graph = readGraphFromFile(inputFilename, representation);
    if (graph == NULL) {
        printf("Erro ao ler o grafo do arquivo.\n");
        return 1;
    }

    writeGraphStatisticsToFile(graph, statsFilename, bfsFilename, dfsFilename, startVertex);

    freeGraph(graph);

    printf("Estatísticas do grafo foram salvas no arquivo %s.\n", statsFilename);
    printf("Resultado da BFS foi salvo no arquivo %s.\n", bfsFilename);
    printf("Resultado da DFS foi salvo no arquivo %s.\n", dfsFilename);

    return 0;
}


