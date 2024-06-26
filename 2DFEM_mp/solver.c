#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "main.h"

int Gauss_elimination(int dim, double *x, double **A, double *b){
    /*ガウス消去法によりAx=b（次元：dim）を解き、解をxに格納する*/
    int i, j, k;
    double* copy_b;
    copy_b = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++) copy_b[i] = b[i];
    /*前進消去*/
    for (i=0; i<dim-1; i++){
        for (j=i+1; j<dim; j++){
            double factor = A[j][i]/A[i][i];
            copy_b[j] -= factor*copy_b[i];
            for  (k=i+1; k<dim; k++) A[j][k] -= factor*A[i][k];
        }
    }
    /*後退代入*/
    copy_b[dim-1] /= A[dim-1][dim-1];
    for (i=dim-2; i>=0; i--){
        for (j=i+1; j<dim; j++) copy_b[i] -= A[i][j] * copy_b[j];
        copy_b[i] /= A[i][i];
    }
    /*解を格納*/
    for (i=0; i<dim; i++){
        x[i] = copy_b[i];
    }
}

int normal_CG(int dim, double *x, double **A, double *b){
    /*共役勾配法によりAx=b（次元：dim）を解き、解をxに格納する*/
    /*次元ステップだけ反復する（Aが正定値対称行列なら収束する）*/
    int i, j;

    double* resid; // 残差ベクトル
    resid = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++){
        resid[i] = b[i];
        for (int j=0; j<dim; j++) resid[i] -= A[i][j]*x[j];
    }
    double* dir; // 探索方向
    dir = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++) dir[i] = resid[i];
    double l2norm_resid; // 残差ベクトルのL2ノルム
    double sum = 0;
    for (int i=0; i<dim; i++){
        sum += resid[i]*resid[i];
    }
    l2norm_resid = sqrt(sum);

    int step = 0; // ステップ数

    double alpha; // 近似解の修正係数
    double rz = 0, pAp = 0; // alphaの因子
    double beta;  // 探索方向の修正係数
    double rz_next = 0; // betaの因子

    double t1, t2;

    t1 = omp_get_wtime();

    double* Adir;
    Adir = (double*)malloc(sizeof(double)*dim); 

    while  (step <= dim){
        /*出力*/
        /*
        if (step%10==0){
            printf("step %4d: L2 norm of the resid %e\n", step, l2norm_resid);
        }
        */

       // スレッド数を2にするとか
        #pragma omp parallel for 
        for (i=0; i<dim; i++){
            Adir[i] = 0;
            #pragma omp parallel for
            for (j=0; j<dim; j++){
                Adir[i] = Adir[i] + A[i][j]*dir[j]; //2重ループにかける
            }
        }

        /*近似解の更新*/
        rz = 0, pAp = 0;
        #pragma omp parallel for reduction(+:pAp)
        for (i=0; i<dim; i++){
            pAp = pAp +  dir[i]*Adir[i];
        }
        #pragma omp parallel for reduction(+:rz)
        for (i=0; i<dim; i++){
            rz = rz + resid[i]*resid[i];
        }
        alpha = rz/pAp;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            x[i] = x[i] + alpha*dir[i];
        }
        
        /*残差ベクトルとそのL2ノルムの更新*/
        #pragma omp parallel for
        for (i=0; i<dim; i++){
                resid[i] = resid[i] - alpha*Adir[i];
        }
        sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for (int i=0; i<dim; i++){
            sum = sum + resid[i]*resid[i];
        }
        l2norm_resid = sqrt(sum);

        /*探索方向の更新*/
        rz_next = 0;
        #pragma omp parallel for reduction(+:rz_next)
        for (i=0; i<dim; i++){
            rz_next = rz_next + resid[i]*resid[i];
        }
        beta = rz_next/rz;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            dir[i] = resid[i] + beta*dir[i];
        }

        /*ステップ数の更新*/
        step++ ;
    }

    t2 = omp_get_wtime();

    int omp_max_threads = omp_get_max_threads();

    printf("solver: CG (not preconditioned)\n");
    printf("using openMP\n");
    printf("CG time = %lf [sec.] \n",t2-t1);
    printf("Problem size = %d (matrix dimension)\n", dim);
    printf("time per size = %lf \n", (t2-t1)/dim);

    free(Adir);
    free(resid);
    free(dir);
}

int diagscaled_CG(int dim, double *x, double **A, double *b){
    /*共役勾配法によりAx=b（次元：dim）を解き、解をxに格納する*/
    /*次元ステップだけ反復する（Aが正定値対称行列なら収束する）*/
    int i, j;

    /*対角スケーリングによる前処理*/
    for (i=0; i<dim; i++){
        double diag = A[i][i];
        for (j=0; j<dim; j++){
            A[i][j] = A[i][j] / diag;
            b[i] = b[i] / diag;
        }
    }

    double* resid; // 残差ベクトル
    resid = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++){
        resid[i] = b[i];
        for (int j=0; j<dim; j++) resid[i] -= A[i][j]*x[j];
    }
    double* dir; // 探索方向
    dir = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++) dir[i] = resid[i];
    double l2norm_resid; // 残差ベクトルのL2ノルム
    double sum = 0;
    for (int i=0; i<dim; i++){
        sum += resid[i]*resid[i];
    }
    l2norm_resid = sqrt(sum);

    int step = 0; // ステップ数

    double alpha; // 近似解の修正係数
    double rz = 0, pAp = 0; // alphaの因子
    double beta;  // 探索方向の修正係数
    double rz_next = 0; // betaの因子

    double t1, t2;

    t1 = omp_get_wtime();

    double* Adir;
    Adir = (double*)malloc(sizeof(double)*dim);

    while  (step <= dim){
        /*出力*/
        /*
        if (step%10==0){
            printf("step %4d: L2 norm of the resid %e\n", step, l2norm_resid);
        }
        */

        for (i=0; i<dim; i++){
            Adir[i] = 0;
            for (j=0; j<dim; j++){
                Adir[i] = Adir[i] + A[i][j]*dir[j];
            }
        }

        /*近似解の更新*/
        rz = 0, pAp = 0;
        #pragma omp parallel for reduction(+:pAp)
        for (i=0; i<dim; i++){
            pAp = pAp +  dir[i]*Adir[i];
        }
        #pragma omp parallel for reduction(+:rz)
        for (i=0; i<dim; i++){
            rz = rz + resid[i]*resid[i];
        }
        alpha = rz/pAp;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            x[i] = x[i] + alpha*dir[i];
        }
        
        /*残差ベクトルとそのL2ノルムの更新*/
        #pragma omp parallel for
        for (i=0; i<dim; i++){
                resid[i] = resid[i] - alpha*Adir[i];
        }
        sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for (int i=0; i<dim; i++){
            sum = sum + resid[i]*resid[i];
        }
        l2norm_resid = sqrt(sum);

        /*探索方向の更新*/
        rz_next = 0;
        #pragma omp parallel for reduction(+:rz_next)
        for (i=0; i<dim; i++){
            rz_next = rz_next + resid[i]*resid[i];
        }
        beta = rz_next/rz;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            dir[i] = resid[i] + beta*dir[i];
        }

        /*ステップ数の更新*/
        step++ ;
    }

    t2 = omp_get_wtime();

    int omp_max_threads = omp_get_max_threads();

    printf("solver: CG (diagonal scaled)\n");
    printf("using openMP\n");
    printf("CG time = %lf [sec.] \n",t2-t1);
    printf("Problem size = %d (matrix dimension)\n", dim);
    printf("time per size = %lf \n", (t2-t1)/dim);

    free(Adir);
    free(resid);
    free(dir);
}

int normal_CG_output(int dim, double *x, double **A, double *b){
    /*共役勾配法によりAx=b（次元：dim）を解き、解をxに格納する*/
    /*次元ステップだけ反復する（Aが正定値対称行列なら収束する）*/
    int i, j;

    double* resid; // 残差ベクトル
    resid = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++){
        resid[i] = b[i];
        for (int j=0; j<dim; j++) resid[i] -= A[i][j]*x[j];
    }
    double* dir; // 探索方向
    dir = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++) dir[i] = resid[i];
    double l2norm_resid; // 残差ベクトルのL2ノルム
    double sum = 0;
    for (int i=0; i<dim; i++){
        sum += resid[i]*resid[i];
    }
    l2norm_resid = sqrt(sum);

    int step = 0; // ステップ数

    double alpha; // 近似解の修正係数
    double rz = 0, pAp = 0; // alphaの因子
    double beta;  // 探索方向の修正係数
    double rz_next = 0; // betaの因子

    double t1, t2;

    t1 = omp_get_wtime();

    double* Adir;
    Adir = (double*)malloc(sizeof(double)*dim); 

    while  (step <= dim){
        /*出力*/
        /*
        if (step%10==0){
            printf("step %4d: L2 norm of the resid %e\n", step, l2norm_resid);
        }
        */

        FILE* file;
        file = fopen("resid_cg.dat", "a");
        fprintf(file, "%e, ", l2norm_resid);
        fclose(file);

        // スレッド数を2にするとか
        #pragma omp parallel for 
        for (i=0; i<dim; i++){
            Adir[i] = 0;
            #pragma omp parallel for
            for (j=0; j<dim; j++){
                Adir[i] = Adir[i] + A[i][j]*dir[j]; //2重ループにかける
            }
        }

        /*近似解の更新*/
        rz = 0, pAp = 0;
        #pragma omp parallel for reduction(+:pAp)
        for (i=0; i<dim; i++){
            pAp = pAp +  dir[i]*Adir[i];
        }
        #pragma omp parallel for reduction(+:rz)
        for (i=0; i<dim; i++){
            rz = rz + resid[i]*resid[i];
        }
        alpha = rz/pAp;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            x[i] = x[i] + alpha*dir[i];
        }
        
        /*残差ベクトルとそのL2ノルムの更新*/
        #pragma omp parallel for
        for (i=0; i<dim; i++){
                resid[i] = resid[i] - alpha*Adir[i];
        }
        sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for (int i=0; i<dim; i++){
            sum = sum + resid[i]*resid[i];
        }
        l2norm_resid = sqrt(sum);

        /*探索方向の更新*/
        rz_next = 0;
        #pragma omp parallel for reduction(+:rz_next)
        for (i=0; i<dim; i++){
            rz_next = rz_next + resid[i]*resid[i];
        }
        beta = rz_next/rz;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            dir[i] = resid[i] + beta*dir[i];
        }

        /*ステップ数の更新*/
        step++ ;
    }

    t2 = omp_get_wtime();

    int omp_max_threads = omp_get_max_threads();

    printf("solver: CG (diagonal scaled)\n");
    printf("using openMP\n");
    printf("CG time = %lf [sec.] \n",t2-t1);
    printf("Problem size = %d (matrix dimension)\n", dim);
    printf("time per size = %lf \n", (t2-t1)/dim);

    FILE* file;
    file = fopen("resid_cg.dat", "a");
    fprintf(file, "%e, ", l2norm_resid);
    fclose(file);

    free(Adir);
    free(resid);
    free(dir);
}

int diagscaled_CG_output(int dim, double *x, double **A, double *b){
    /*共役勾配法によりAx=b（次元：dim）を解き、解をxに格納する*/
    /*次元ステップだけ反復する（Aが正定値対称行列なら収束する）*/
    int i, j;

    /*対角スケーリングによる前処理*/
    for (i=0; i<dim; i++){
        double diag = A[i][i];
        for (j=0; j<dim; j++){
            A[i][j] = A[i][j] / diag;
            b[i] = b[i] / diag;
        }
    }

    double* resid; // 残差ベクトル
    resid = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++){
        resid[i] = b[i];
        for (int j=0; j<dim; j++) resid[i] -= A[i][j]*x[j];
    }
    double* dir; // 探索方向
    dir = (double*)malloc(sizeof(double)*dim);
    for (i=0; i<dim; i++) dir[i] = resid[i];
    double l2norm_resid; // 残差ベクトルのL2ノルム
    double sum = 0;
    for (int i=0; i<dim; i++){
        sum += resid[i]*resid[i];
    }
    l2norm_resid = sqrt(sum);

    int step = 0; // ステップ数

    double alpha; // 近似解の修正係数
    double rz = 0, pAp = 0; // alphaの因子
    double beta;  // 探索方向の修正係数
    double rz_next = 0; // betaの因子

    double t1, t2;

    t1 = omp_get_wtime();

    double* Adir;
    Adir = (double*)malloc(sizeof(double)*dim);

    while  (step <= dim){
        /*出力*/
        /*
        if (step%10==0){
            printf("step %4d: L2 norm of the resid %e\n", step, l2norm_resid);
        }
        */

        FILE* file;
        file = fopen("resid_dscg.dat", "a");
        fprintf(file, "%e, ", l2norm_resid);
        fclose(file);

        for (i=0; i<dim; i++){
            Adir[i] = 0;
            for (j=0; j<dim; j++){
                Adir[i] = Adir[i] + A[i][j]*dir[j];
            }
        }

        /*近似解の更新*/
        rz = 0, pAp = 0;
        #pragma omp parallel for reduction(+:pAp)
        for (i=0; i<dim; i++){
            pAp = pAp +  dir[i]*Adir[i];
        }
        #pragma omp parallel for reduction(+:rz)
        for (i=0; i<dim; i++){
            rz = rz + resid[i]*resid[i];
        }
        alpha = rz/pAp;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            x[i] = x[i] + alpha*dir[i];
        }
        
        /*残差ベクトルとそのL2ノルムの更新*/
        #pragma omp parallel for
        for (i=0; i<dim; i++){
                resid[i] = resid[i] - alpha*Adir[i];
        }
        sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for (int i=0; i<dim; i++){
            sum = sum + resid[i]*resid[i];
        }
        l2norm_resid = sqrt(sum);

        /*探索方向の更新*/
        rz_next = 0;
        #pragma omp parallel for reduction(+:rz_next)
        for (i=0; i<dim; i++){
            rz_next = rz_next + resid[i]*resid[i];
        }
        beta = rz_next/rz;
        #pragma omp parallel for
        for (i=0; i<dim; i++){
            dir[i] = resid[i] + beta*dir[i];
        }

        /*ステップ数の更新*/
        step++ ;
    }

    t2 = omp_get_wtime();

    int omp_max_threads = omp_get_max_threads();

    printf("solver: CG (diagonal scaled)\n");
    printf("using openMP\n");
    printf("CG time = %lf [sec.] \n",t2-t1);
    printf("Problem size = %d (matrix dimension)\n", dim);
    printf("time per size = %lf \n", (t2-t1)/dim);

    FILE* file;
    file = fopen("resid_dscg.dat", "a");
    fprintf(file, "%e, ", l2norm_resid);
    fclose(file);


    free(Adir);
    free(resid);
    free(dir);
}