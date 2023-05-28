#ifndef _cond_cont_h_
#define _cond_cont_h_
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                    /** CONDICIONES DE CONTORNO **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

void cond_cont ( double **K, double *F, struct adyacencias *tabla_ady,
                 double **coord, int nodos, int tipo_elem, int nodos_elem,
                 int FP , int  *f_punt_ind, double **f_punt_val,
                 int FPP, int **f_pre_ind , double  *f_pre_val ,
                 int FS , int **f_sup_ind , double **f_sup_val ,
                 int DL , int  *d_lin_ind , double **d_lin_val ,
                 int DS , int **d_sup_ind , double **d_sup_val , int dim);

void dirichlet ( double **K, double *F, struct adyacencias *tabla_ady,
                 int DL, int DS, int *d_lin_ind, double **d_lin_val,
                 int **d_sup_ind, double **d_sup_val, int nodos,
                 int tipo_elem, int dim );

void neumann ( double **coord, double *F, int FP, int FPP, int FS,
               int *f_punt_ind, double **f_punt_val, int **f_pre_ind,
               double *f_pre_val, int **f_sup_ind, double **f_sup_val,
               int tipo_elem, int nodos_elem, int dim );

void dirich_DL ( double **K, double *F, struct adyacencias *tabla_ady,
                 int *d_lin_ind, double **d_lin_val, int DL, int nodos,
                 int dim );

void dirich_DS ( double **K, double *F, struct adyacencias *tabla_ady,
                 int **d_sup_ind, double **d_sup_val, int DS, int nodos,
                 int tipo_elem, int dim );

void neu_FS  ( double *F, int **f_sup_ind, double **f_sup_val, int FS,
               double **coord, int tipo_elem, int nodos_elem, int dim );

void neu_FP  ( double *F, int *f_punt_ind, double **f_punt_val, int FP,
               int dim );

void neu_FPP ( double *F, int **f_pre_ind, double *f_pre_val, int FPP,
               double **coord, int tipo_elem, int nodos_elem, int dim );

void FPP_3 ( double *F, int **f_pre_ind, double *f_pre_val, int FPP,
             double **coord, int dim );

void FPP_4 ( double *F, int **f_pre_ind, double *f_pre_val, int FPP,
             double **coord, int dim );

void FS_3  ( double *F, int **f_sup_ind, double **f_sup_val, int FS,
             double **coord, int dim );

void FS_4  ( double *F, int **f_sup_ind, double **f_sup_val, int FS,
             double **coord, int dim );

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

