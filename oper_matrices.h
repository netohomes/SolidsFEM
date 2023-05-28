#ifndef _oper_matrices_h_
#define _oper_matrices_h_
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                     /**  OPERACIONES MATRICIALES Y VECTORIALES  **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

void multiplica ( double **C, double **A, double **B, int rengA, int col, int colB );
void mult_mat_vector ( double *C, double **A, double *B, int rengA, int colA, int rengB );
void transponer_mat ( double **A, double **B, int reng, int col );
void mult_mat_escalar ( double **A, double b, int reng, int col);
void mult_vec_escalar ( double *v, double b, int reng );
void suma_mat ( double **A, double **B, int reng, int col);
void suma_vec ( double *a, double *b, int reng );

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
