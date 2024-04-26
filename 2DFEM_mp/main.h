#ifndef _MAIN_H_
#define _MAIN_H_

#define L 50 // 一辺の分割個数 ここを入力！！
#define NUMOFELEMS L*L*2 // 有限要素の数
#define NUMOFNODES (L+1)*(L+1) // 節点の数
#define N L*4 /*境界節点の数*/
#define M (L-1)*(L-1) /*非境界節点の数*/

int input_connectivity(int **cnn);
int input_node_coord(double *node_coord_x, double *node_coord_y);
int input_border(double *u, double *node_coord_x, double *node_coord_y);
int input_nodedata(int *border_nodes, int *nonborder_nodes);

int construct_matrix(double **A, double *f, int **cnn, double *node_coord_x, double *node_coord_y, double F, int *border_nodes, int *nonborder_nodes, double *u);

int Gauss_elimination(int dim, double *x, double **A, double *b);
int normal_CG(int dim, double *x, double **A, double *b);
int diagscaled_CG(int dim, double *x, double **A, double *b);
int normal_CG_output(int dim, double *x, double **A, double *b);
int diagscaled_CG_output(int dim, double *x, double **A, double *b);

#endif // _MAIN_H_