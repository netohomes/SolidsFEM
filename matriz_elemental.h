#ifndef _matriz_elemental_h_
#define _matriz_elemental_h_
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                                /** MATRICES ELEMENTALES **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

void rigid_elem ( double **K_elem, double *f_elem, double **D, double **coord, double **B,
                  double **BT, double *det, double **pg, double **wpg, double **jac,
                  double **derin, double **Baux, double **K_elemaux, int dim, int tipo_elem,
                  int nodos_elem, int pg_elem, double **NT, double *g, double *f_elem_aux,
                  int m, double E, double v );

void fuer_elem ( double *f_elem_aux, double **NT, double *b, double p, double n, double c, int dim,
                 int tipo_elem, int nodos_elem );

void mat_N ( double **NT, double p, double n, double c, int dim, int tipo_elem, int nodos_elem,
            int rNT, int cNT );

void mat_B ( double **B, int dim, int tipo_elem, int nodos_elem, double **inv_jac, double p,
             double n, double c );

void tam_mat ( int tipo_elem, int dim, int nodos_elem, int *rB, int *cB, int *rD, int *cD, int *rf );

void sum_B ( double **B, double *v2, int dim, int j );

void mat_jacobiana ( double **jacobiano, double **coord, double p, double n, double c, int nodos_elem,
                     int tipo_elem, int dim, double **derin );

void deriv_nat ( double **derin, double p, double n, double c, int tipo_elem );

double determinante ( double **A, int dim );

void inversa ( double **A, double **inversa_A, int n, double det );

void deriv_Ni( double *vect_aux, double p, double n, double c, int tipo_elem, int dim, int j );

void PG_k ( double **pg, double **wpg, int dim, int k, double *p, double *n, double *c,
            double *wp,  double *wn,   double *wc );

void constitutiva ( double **A, double E, double v, int dim );

double N ( int tipo_elem, int i, int dim, double p, double n, double c );

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

