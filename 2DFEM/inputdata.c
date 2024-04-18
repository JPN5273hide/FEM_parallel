#include <stdio.h>
#include <stdlib.h>
#include "main.h"

int input_connectivity(int **cnn){
    /*コネクティビティを設定*/
    int **array_cnn;
    array_cnn = (int**)malloc(sizeof(int*)*NUMOFELEMS);
    for (int i=0; i<NUMOFELEMS; i++) *(array_cnn+i) = (int*)malloc(sizeof(int)*3);
    for (int elem=0; elem<NUMOFELEMS; elem++){
        int blocknum = elem/2;
        int i = blocknum/L;
        int j = blocknum%L;
        if (elem%2==0){
            array_cnn[elem][0] = i*(L+1)+j;
            array_cnn[elem][1] = (i+1)*(L+1)+j+1;
            array_cnn[elem][2] = (i+1)*(L+1)+j;
        }
        else if (elem%2==1){
            array_cnn[elem][0] = i*(L+1)+j;
            array_cnn[elem][1] = i*(L+1)+j+1;
            array_cnn[elem][2] = (i+1)*(L+1)+j+1;
        }
    }
    for (int elem=0; elem<NUMOFELEMS; elem++){
        for (int elemnode=0; elemnode<3; elemnode++){
            cnn[elem][elemnode] = array_cnn[elem][elemnode];
        }
    }

    free(array_cnn);

    return 0;
}

int input_node_coord(double *node_coord_x, double *node_coord_y){
    /*節点の座標を設定*/
    for (int node=0; node<NUMOFNODES; node++){
        *(node_coord_x+node) = (node%(L+1))/(double)L;
        *(node_coord_y+node) = (node/(L+1))/(double)L;
    }
    return 0;
}

int input_border(double *u, double *node_coord_x, double *node_coord_y){
    /*境界条件を設定*/
    for (int i=0; i<NUMOFNODES; i++){
        double x = *(node_coord_x+i);
        double y = *(node_coord_y+i);
        if (x==0) *(u+i) = y*y;
        else if (y==0) *(u+i) = x*x;
        else if (x==1) *(u+i) = 1+y*y;
        else if (y==1) *(u+i) = 1+x*x;
        else *(u+i) = 0;
    }

    return 0;
}

int input_nodedata(int *border_nodes, int *nonborder_nodes){
    /*境界節点と非境界節点を識別*/
    int iter_to_N = 0, iter_to_M = 0;
    int node = 0;
    for(node=0; node<NUMOFNODES; node++){
        int quo = node/(int)(L+1);
        int rem = node%(int)(L+1);
        if (quo == 0){
            border_nodes[iter_to_N] = node;
            iter_to_N++;
        }
        else if (quo == L){
            border_nodes[iter_to_N] = node;
            iter_to_N++;
        }
        else if (rem == 0 ){
            border_nodes[iter_to_N] = node;
            iter_to_N++;
        }
        else if (rem == L){
            border_nodes[iter_to_N] = node;
            iter_to_N++;
        }
        else {
            nonborder_nodes[iter_to_M] = node;
            iter_to_M++;
        }
    }

    return 0;
}