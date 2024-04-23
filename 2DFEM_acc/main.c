#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "main.h"

int main(void){

    int** cnn;
    cnn = (int**)malloc(sizeof(int*)*NUMOFELEMS);
    for (int i=0; i<NUMOFELEMS; i++) *(cnn+i) = (int*)malloc(sizeof(int)*3);
    double* node_coord_x;
    node_coord_x = (double*)malloc(sizeof(double)*NUMOFNODES);
    double* node_coord_y;
    node_coord_y = (double*)malloc(sizeof(double)*NUMOFNODES);
    input_connectivity(cnn);
    input_node_coord(node_coord_x, node_coord_y);

    double* u;
    u = (double*)malloc(sizeof(double)*NUMOFNODES);
    input_border(u, node_coord_x, node_coord_y);

    int* border_nodes;
    border_nodes = (int*)malloc(sizeof(int)*N);
    int* non_bordernodes;
    non_bordernodes = (int*)malloc(sizeof(int)*M);
    input_nodedata(border_nodes, non_bordernodes);

    double** A;
    A = (double**)malloc(sizeof(double*)*M);
    for (int i=0; i<M; i++) *(A+i) = (double*)malloc(sizeof(double)*M);
    double* f;
    f = (double*)malloc(sizeof(double)*M);
    construct_matrix(A, f, cnn, node_coord_x, node_coord_y, -4.0, border_nodes, non_bordernodes, u);

    free(cnn);
    free(border_nodes);
    free(non_bordernodes);

    double x[M];
    for (int i=0; i<M; i++) x[i] = 1.0;
    normal_CG(M, x, A, f);

    free(node_coord_x);
    free(node_coord_y);
    free(u);
    free(A);
    free(f);
    return 0;
}