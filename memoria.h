#ifndef _memoria_h_
#define _memoria_h_
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                                    /**  ASIGNA MEMORIA  **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

double **mem_matriz ( int n, int  m, char *nombre );
int **mem_matriz_int ( int n, int  m, char *nombre );
double *mem_vector ( int n, char  *nombre );
int *mem_vector_int ( int n, char  *nombre );
struct elemental *mem_vector_elem ( int n, char  *nombre );
void redim_vector_elem ( struct elemental *a, int n, char  *nombre );
struct adyacencias *mem_vector_ady ( int n, char  *nombre );
void redim_vector_ady ( struct adyacencias *A, int m,int n, char *nombre );
void redimensiona_vector ( int **A, int m, int n, char *nombre );
void redim_matriz ( double **A, int n, int  m, char *nombre );
void inicializa_matriz ( double **A, int reng, int col );
void inicializa_vector ( double *A, int reng );
void inicializa_vector_int ( int *A, int reng );
void inicializa_mat_identidad ( double **A, int n );
void libera_matriz ( double **a, int reng );
void libera_matriz_int ( int **a, int reng );
void libera_vector ( double *v );
void libera_vector_int ( int *v );
void libera_elem ( struct elemental *a, int n_elem );
//void libera_tipo_2 ( struct tipo_2 *a, int nodos );
void libera_mem1 ( double *T, double **q_pg, double **q_pg2, double **q_nodos, double **error_PG,
                   int dimension, double *det );
void libera_mem2 ( void );

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
