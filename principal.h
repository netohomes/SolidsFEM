#ifndef _principal_h_
#define _principal_h_
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

        /** MECÁNICA DE SÓLIDOS RESUELTO POR ELEMENTOS FINITOS **/
                             // Para 3D

    // Estructuras
    struct elemental{
        int *conec;
        double E;
        double v;
    };
    struct adyacencias{
        int *ind;
        int cols;
    };

    // Variables
    int dim, tipo_elem, pg_elem, n_elem, nodos, nodos_elem, itera;
    struct elemental   *elemento;
    struct adyacencias *tabla_ady;
    double **coord, *det, **K, *F, *U, **ESF, error;
    int FP, FPP, FS, FC, DL, DS;
    int     *f_punt_ind, **f_pre_ind, **f_sup_ind,  *f_cpo_ind;
    double **f_punt_val,  *f_pre_val, **f_sup_val, **f_cpo_val;
    int      *d_lin_ind, **d_sup_ind;
    double  **d_lin_val, **d_sup_val;

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}


