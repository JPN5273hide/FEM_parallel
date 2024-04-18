#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "main.h"

int construct_matrix(double **A, double *f, int **cnn, double *node_coord_x, double *node_coord_y, double F, int *border_nodes, int *nonborder_nodes, double *u){
    /*拡大係数行列と拡大定数項ベクトルを初期化*/
    double** largeA;
    largeA = (double**)malloc(sizeof(double*)*NUMOFNODES);
    for (int i=0; i<NUMOFNODES; i++) *(largeA+i) = (double*)malloc(sizeof(double)*NUMOFNODES);
    for (int i=0; i<NUMOFNODES; i++) for (int j=0; j<NUMOFNODES; j++) largeA[i][j] = 0;
    double* largef;
    largef = (double*)malloc(sizeof(double)*NUMOFNODES);
    for (int i=0; i<NUMOFNODES; i++) largef[i] = 0;
    
    /*拡大係数行列, 拡大定数項ベクトルの構成*/
    for (int elem=0; elem<NUMOFELEMS; elem++){
        /*座標データの読み込み*/
        int node[3];
        for (int i=0; i<3; i++) node[i] = cnn[elem][i];
        double x0 = node_coord_x[node[0]];
        double x1 = node_coord_x[node[1]];
        double x2 = node_coord_x[node[2]];
        double y0 = node_coord_y[node[0]];
        double y1 = node_coord_y[node[1]];
        double y2 = node_coord_y[node[2]];

        /*要素面積*/
        double D = x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1);
        double S = fabs(D)/2;

        /*要素係数行列の構成*/
        double elemA[3][3];
        double b[3] = {y1-y2, y2-y0, y0-y1};
        double c[3] = {x2-x1, x0-x2, x1-x0};
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) elemA[i][j] = (b[j]*b[i]+c[j]*c[i])/(4.0*S);

        /*要素定数項ベクトルの構成*/
        double elemf[3];
        for (int i=0; i<3; i++) elemf[i] = F*S/3.0;

        /*拡大係数行列, 拡大定数項ベクトルに組み込む*/
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) largeA[node[i]][node[j]] += elemA[i][j];
        for (int i=0; i<3; i++) largef[node[i]] += elemf[i];
    }

    /*境界条件の考慮*/  
    for (int i=0; i<M; i++) for (int j=0; j<M; j++) A[i][j] = largeA[nonborder_nodes[i]][nonborder_nodes[j]];
    for (int i=0; i<M; i++) f[i] = largef[nonborder_nodes[i]];
    for (int i=0; i<N; i++){
        for (int j=0; j<M; j++){
            f[j] -= largeA[nonborder_nodes[j]][border_nodes[i]]*u[border_nodes[i]];
        }
    } 

    free(largeA);
    free(largef);

    return 0;
}