#ifndef _matriz_global_h_
#define _matriz_global_h_
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                    /** MATRIZ DE RIGIDEZ GLOBAL **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

    void ady ( struct elemental *elem, int n_elem, int nodos,
               int nodos_elem );

    void rigid ( struct adyacencias *tabla_ady, struct elemental *elem,
                 int *f_cpo_ind, double **f_cpo_val, int tipo_elem,
                 int nodos_elem, int pg_elem, int dim );

    void P_G ( double **pg, double **wpg, int tipo_elem, int dim );

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
