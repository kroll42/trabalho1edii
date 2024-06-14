// grafo.h
#define GRAFO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LIST 0
#define MATRIX 1

typedef struct AdjListNode {
    int dest;
    struct AdjListNode* next;
} AdjListNode;

typedef struct AdjList {
    AdjListNode *head;
} AdjList;

typedef struct Graph {
    int numVertices;
    int numEdges;
    AdjList* array;
    int** adjMatrix; // Matriz de adjacência
    int representation; // 0 para Lista, 1 para Matriz
} Graph;


Graph* createGraph(int numVertices, int representation) {
    Graph* graph = (Graph*) malloc(sizeof(Graph));
    graph->numVertices = numVertices;
    graph->numEdges = 0;
    graph->representation = representation;

    if (representation == LIST) {
        graph->array = (AdjList*) malloc(numVertices * sizeof(AdjList));
        for (int i = 0; i < numVertices; ++i) {
            graph->array[i].head = NULL;
        }
    } else if (representation == MATRIX) {
        graph->adjMatrix = (int**) malloc(numVertices * sizeof(int*));
        for (int i = 0; i < numVertices; i++) {
            graph->adjMatrix[i] = (int*) calloc(numVertices, sizeof(int));
        }
    }

    return graph;
}

AdjListNode* createNode(int dest) {
    AdjListNode* newNode = (AdjListNode*) malloc(sizeof(AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}

void addEdge(Graph* graph, int src, int dest) {
    if (graph->representation == LIST) {
        AdjListNode* newNode = createNode(dest);
        newNode->next = graph->array[src].head;
        graph->array[src].head = newNode;

        newNode = createNode(src);
        newNode->next = graph->array[dest].head;
        graph->array[dest].head = newNode;
    } else if (graph->representation == MATRIX) {
        graph->adjMatrix[src][dest] = 1;
        graph->adjMatrix[dest][src] = 1;
    }

    graph->numEdges++;
}

Graph* readGraphFromFile(const char* filename, int representation) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Erro ao abrir o arquivo!\n");
        return NULL;
    }

    int numVertices;
    fscanf(file, "%d", &numVertices);

    Graph* graph = createGraph(numVertices, representation);

    int src, dest;
    while (fscanf(file, "%d %d", &src, &dest) != EOF) {
        addEdge(graph, src, dest);
    }

    fclose(file);
    return graph;
}
void DFS(Graph* graph, int startVertex, int* parent, int* level, FILE* file) {
    int* visited = (int*) calloc(graph->numVertices, sizeof(int));
    if (visited == NULL) {
        printf("Erro de alocação de memória para visited\n");
        return;
    }

    int* stack = (int*) malloc(graph->numVertices * sizeof(int));
    if (stack == NULL) {
        printf("Erro de alocação de memória para stack\n");
        free(visited);
        return;
    }

    int top = -1;
    stack[++top] = startVertex;
    parent[startVertex] = -1; // Raiz não tem pai
    level[startVertex] = 0;

    while (top >= 0) {
        int v = stack[top--];

        if (!visited[v]) {
            visited[v] = 1;
            fprintf(file, "Vértice %d: Pai: %d, Nível: %d\n", v, parent[v], level[v]);

            // Obter adjacentes
            if (graph->representation == LIST) {
                AdjListNode* adj = graph->array[v].head;
                while (adj) {
                    if (!visited[adj->dest]) {
                        stack[++top] = adj->dest;
                        parent[adj->dest] = v;
                        level[adj->dest] = level[v] + 1;
                    }
                    adj = adj->next;
                }
            } else if (graph->representation == MATRIX) {
                for (int i = 0; i < graph->numVertices; i++) {
                    if (graph->adjMatrix[v][i] && !visited[i]) {
                        stack[++top] = i;
                        parent[i] = v;
                        level[i] = level[v] + 1;
                    }
                }
            }
        }
    }

    free(visited);
    free(stack);
}

void BFS(Graph* graph, int startVertex, int* parent, int* level, FILE* file) {
    int* visited = (int*) calloc(graph->numVertices, sizeof(int));
    if (visited == NULL) {
        printf("Erro de alocação de memória para visited\n");
        return;
    }

    int* queue = (int*) malloc(graph->numVertices * sizeof(int));
    if (queue == NULL) {
        printf("Erro de alocação de memória para queue\n");
        free(visited);
        return;
    }

    int front = 0, rear = 0;
    queue[rear++] = startVertex;
    parent[startVertex] = -1; // Raiz não tem pai
    level[startVertex] = 0;
    visited[startVertex] = 1;

    while (front < rear) {
        int v = queue[front++];

        fprintf(file, "Vértice %d: Pai: %d, Nível: %d\n", v, parent[v], level[v]);

        // Obter adjacentes
        if (graph->representation == LIST) {
            AdjListNode* adj = graph->array[v].head;
            while (adj) {
                if (!visited[adj->dest]) {
                    queue[rear++] = adj->dest;
                    parent[adj->dest] = v;
                    level[adj->dest] = level[v] + 1;
                    visited[adj->dest] = 1;
                }
                adj = adj->next;
            }
        } else if (graph->representation == MATRIX) {
            for (int i = 0; i < graph->numVertices; i++) {
                if (graph->adjMatrix[v][i] && !visited[i]) {
                    queue[rear++] = i;
                    parent[i] = v;
                    level[i] = level[v] + 1;
                    visited[i] = 1;
                }
            }
        }
    }

    free(visited);
    free(queue);
}
// Função de comparação para qsort
int cmpfunc(const void* a, const void* b) {
    return (*(int*)b - *(int*)a);
}

void calculateDegreeStatistics(Graph* graph, int* minDegree, int* maxDegree, double* avgDegree, int* medianDegree) {
    int* degrees = (int*) calloc(graph->numVertices, sizeof(int));
    int totalDegree = 0;

    if (graph->representation == LIST) {
        for (int i = 0; i < graph->numVertices; i++) {
            AdjListNode* node = graph->array[i].head;
            while (node) {
                degrees[i]++;
                node = node->next;
            }
            totalDegree += degrees[i];
        }
    } else if (graph->representation == MATRIX) {
        for (int i = 0; i < graph->numVertices; i++) {
            for (int j = 0; j < graph->numVertices; j++) {
                if (graph->adjMatrix[i][j]) {
                    degrees[i]++;
                }
            }
            totalDegree += degrees[i];
        }
    }

    // Calcular grau mínimo, máximo e médio
    *minDegree = degrees[0];
    *maxDegree = degrees[0];
    for (int i = 0; i < graph->numVertices; i++) {
        if (degrees[i] < *minDegree) *minDegree = degrees[i];
        if (degrees[i] > *maxDegree) *maxDegree = degrees[i];
    }
    *avgDegree = (double)totalDegree / graph->numVertices;

    // Calcular a mediana do grau
    qsort(degrees, graph->numVertices, sizeof(int), cmpfunc);
    if (graph->numVertices % 2 == 0)
        *medianDegree = (degrees[graph->numVertices / 2 - 1] + degrees[graph->numVertices / 2]) / 2;
    else
        *medianDegree = degrees[graph->numVertices / 2];

    free(degrees);
}

void DFSUtil(Graph* graph, int v, int visited[], int* component, int* size, int componentIndex, FILE* file) {
    visited[v] = 1;
    component[(*size)++] = v;

    AdjListNode* node = graph->array[v].head;
    while (node) {
        if (!visited[node->dest]) {
            DFSUtil(graph, node->dest, visited, component, size, componentIndex, file);
        }
        node = node->next;
    }
}

void findConnectedComponents(Graph* graph, FILE* file) {
    int* visited = (int*) calloc(graph->numVertices, sizeof(int));
    int* component = (int*) malloc(graph->numVertices * sizeof(int));
    int size;
    int* componentSizes = (int*) malloc(graph->numVertices * sizeof(int));
    int numComponents = 0;

    for (int v = 0; v < graph->numVertices; v++) {
        if (!visited[v]) {
            size = 0;
            DFSUtil(graph, v, visited, component, &size, numComponents, file);
            componentSizes[numComponents++] = size;

            fprintf(file, "Componente %d (Tamanho %d):\n", numComponents, size);
            for (int i = 0; i < size; i++) {
                fprintf(file, "%d ", component[i]);
            }
            fprintf(file, "\n");
        }
    }

    // Ordenar componentes por tamanho decrescente
    qsort(componentSizes, numComponents, sizeof(int), cmpfunc);
    fprintf(file, "Numero de componentes conexos: %d\n", numComponents);
    for (int i = 0; i < numComponents; i++) {
        fprintf(file, "Componente %d: %d vertices\n", i + 1, componentSizes[i]);
    }

    free(visited);
    free(component);
    free(componentSizes);
}

void writeGraphStatisticsToFile(Graph* graph, const char* statsFilename, const char* bfsFilename, const char* dfsFilename, int startVertex) {
    FILE *statsFile = fopen(statsFilename, "w");
    if (statsFile == NULL) {
        printf("Erro ao abrir o arquivo de estatisticas para escrita!\n");
        return;
    }

    int minDegree, maxDegree, medianDegree;
    double avgDegree;
    calculateDegreeStatistics(graph, &minDegree, &maxDegree, &avgDegree, &medianDegree);

    fprintf(statsFile, "Estatisticas do Grafo:\n");
    fprintf(statsFile, "Numero de vertices: %d\n", graph->numVertices);
    fprintf(statsFile, "Numero de arestas: %d\n", graph->numEdges);
    fprintf(statsFile, "Grau minimo: %d\n", minDegree);
    fprintf(statsFile, "Grau maximo: %d\n", maxDegree);
    fprintf(statsFile, "Grau medio: %.2f\n", avgDegree);
    fprintf(statsFile, "Mediana de grau: %d\n\n", medianDegree);

    fclose(statsFile);

    // Realizar DFS
    FILE *dfsFile = fopen(dfsFilename, "w");
    if (dfsFile == NULL) {
        printf("Erro ao abrir o arquivo DFS para escrita!\n");
        return;
    }
    fprintf(dfsFile, "DFS a partir do vertice %d:\n", startVertex);
    int* parentDFS = (int*) malloc(graph->numVertices * sizeof(int));
    int* levelDFS = (int*) malloc(graph->numVertices * sizeof(int));
    if (parentDFS == NULL || levelDFS == NULL) {
        printf("Erro de alocação de memória para DFS!\n");
        fclose(dfsFile);
        return;
    }
    memset(parentDFS, -1, graph->numVertices * sizeof(int));
    memset(levelDFS, -1, graph->numVertices * sizeof(int));
    DFS(graph, startVertex, parentDFS, levelDFS, dfsFile);
    fprintf(dfsFile, "\n");
    fclose(dfsFile);

    free(parentDFS);
    free(levelDFS);

    // Realizar BFS
    FILE *bfsFile = fopen(bfsFilename, "w");
    if (bfsFile == NULL) {
        printf("Erro ao abrir o arquivo BFS para escrita!\n");
        return;
    }
    fprintf(bfsFile, "BFS a partir do vertice %d:\n", startVertex);
    int* parentBFS = (int*) malloc(graph->numVertices * sizeof(int));
    int* levelBFS = (int*) malloc(graph->numVertices * sizeof(int));
    if (parentBFS == NULL || levelBFS == NULL) {
        printf("Erro de alocação de memória para BFS!\n");
        fclose(bfsFile);
        return;
    }
    memset(parentBFS, -1, graph->numVertices * sizeof(int));
    memset(levelBFS, -1, graph->numVertices * sizeof(int));
    BFS(graph, startVertex, parentBFS, levelBFS, bfsFile);
    fprintf(bfsFile, "\n");
    fclose(bfsFile);

    free(parentBFS);
    free(levelBFS);
}


void freeGraph(Graph* graph) {
    if (graph->representation == LIST) {
        for (int i = 0; i < graph->numVertices; i++) {
            AdjListNode* current = graph->array[i].head;
            while (current) {
                AdjListNode* next = current->next;
                free(current);
                current = next;
            }
        }
        free(graph->array);
    } else if (graph->representation == MATRIX) {
        for (int i = 0; i < graph->numVertices; i++) {
            free(graph->adjMatrix[i]);
        }
        free(graph->adjMatrix);
    }
    free(graph);
}

void printAdjMatrix(Graph* graph) {
    if (graph->representation != MATRIX) {
        printf("Grafo nao esta representado como matriz de adjacencia.\n");
        return;
    }

    printf("Matriz de Adjacencia:\n");
    for (int i = 0; i < graph->numVertices; i++) {
        for (int j = 0; j < graph->numVertices; j++) {
            printf("%d ", graph->adjMatrix[i][j]);
        }
        printf("\n");
    }
}


