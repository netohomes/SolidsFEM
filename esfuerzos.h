#ifndef _esfuerzos_h_
#define _esfuerzos_h_
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                  /** CALCULA LOS ESFUERZOS **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

void esfuerzos ( double *U, double **coord, struct elemental *elem, int dim,
                 int nodos, int n_elem, int tipo_elem, int nodos_elem );

void def ( double *e, double **coord, double *u, double **B, double **jac,
           double **inv_jac, double **derin, int tipo_elem, int nodos_elem,
           int rB, int cB, int dim, double p, double n, double c );

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
