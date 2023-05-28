//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                  /** CALCULA LOS ESFUERZOS **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "principal.h"
#include "esfuerzos.h"
#include "memoria.h"
#include "matriz_global.h"
#include "matriz_elemental.h"
#include "oper_matrices.h"
#include <stdio.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void esfuerzos ( double *U, double **coord, struct elemental *elem, int dim,
                 int nodos, int n_elem, int tipo_elem, int nodos_elem )
{
    /** Descripción **/
    // Calcula los esfuerzos en los puntos de Gauss

    /** Variables de entrada **/
    // {U}          Vector de desplazamientos
    // n_elem       Número de elementos
    // nodos_elem   Nodos del elemento
    // dim          Dimensión
    // [ESF]        Vectores de esfuerzos en los P.G. para cada elemento.

    /** Variables locales **/
    // e        Vector de deformaciones
    // p, n, c  Coordendas P.G. en el espacio normalizado

    int rB, cB, rD, cD, rf, a, m, i, j, k, b;
    double p, n, c, wp, wn, wc, **pg, **wpg, **B, **jac, **inv_jac, *u, *e;
    double **D, **derin, **coord_aux, *esf, E, v;

    tam_mat ( tipo_elem, dim, nodos_elem, &rB, &cB, &rD, &cD, &rf );

    a       = nodos_elem;
    ESF     = mem_matriz ( pg_elem * n_elem, rD, "esfuerzos." );
    u       = mem_vector ( cB, "despl. de un elemento." );
    B       = mem_matriz ( rB, cB, "B." );
    jac     = mem_matriz ( dim, dim, "del Jacobiano." );
    inv_jac = mem_matriz ( dim, dim,  "del inv_jacobiano." );
    e       = mem_vector ( rD, "deformaciones." );
    esf     = mem_vector ( rD, "esfuerzos de un elem." );
    D       = mem_matriz ( rD, cD, "Constitutiva." );
    pg      = mem_matriz ( pg_elem, dim, "puntos de Gauss." );
    wpg     = mem_matriz ( pg_elem, dim, "pesos de P.G." );
    derin   = mem_matriz ( dim, a,  "derivadas naturales." );
    coord_aux = mem_matriz ( a, dim, "coord_aux para elems.");

    P_G ( pg, wpg, tipo_elem, dim );

    for ( m = 0; m < n_elem; m++ ){

        E = elem[ m ].E;
        v = elem[ m ].v;

        constitutiva ( D, E, v, dim );

        for ( i = 0; i < nodos_elem; i++ ){
            b = i * dim;
            for ( j = 0; j < dim; j++ ){
                a = elem[ m ].conec[ i ];
                coord_aux[ i ][ j ] = coord[ a ][ j ];
                a *= dim;
                u[ b + j ] = U[ a + j ];
            }
        }

        //printf ( "\n\nu" );
        //for ( i = 0; i < cB; i++ ) printf ( "\n%lf", u[ i ] );

        for ( k = 0; k < pg_elem; k++ ){

            PG_k ( pg, wpg, dim, k, &p, &n, &c, &wp, &wn, &wc );

            def ( e, coord_aux, u, B, jac, inv_jac, derin, tipo_elem,
                  nodos_elem, rB, cB, dim, p, n, c );

            /*printf ( "\n\nDef." );
            for ( i = 0; i < rD; i++ ) printf ( "\n%lf", e[ i ] );*/

            inicializa_vector ( esf, rD );

            mult_mat_vector ( esf, D, e, rD, cD, cD );

            /*printf ( "\n\nEsf." );
            for ( i = 0; i < rD; i++ ) printf ( "\n%lf", esf[ i ] );*/

            for ( i = 0; i < rD; i++ ){

                ESF[ m * pg_elem + k ][ i ] = esf [ i ];

            }

        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void def ( double *e, double **coord, double *u, double **B, double **jac,
           double **inv_jac, double **derin, int tipo_elem, int nodos_elem,
           int rB, int cB, int dim, double p, double n, double c )
{
    /** Descripción **/
    // Calcula el vector de deformaciones para un P.G

    /** Variables de entrada **/
    // {U}      Vector de desplazamientos
    // [B]      Matriz para guardar a B.
    // e        Vector de deformaciones

    double det_jac;
    //int a, b;

    inicializa_matriz ( inv_jac, dim, dim );
    inicializa_matriz ( B, rB, cB );
    inicializa_vector ( e, rB );

    mat_jacobiana ( jac, coord, p, n, c, nodos_elem, tipo_elem, dim, derin );

        /*printf ( "\njacobiano\n");
        for ( a = 0; a < dim; a++ ){
            for ( b = 0; b < dim; b++ ){
                printf ( "%lf\t", jac[a][b] );
                if ( b == dim - 1 ) printf ( "\n" );
            }
        }*/

    det_jac = determinante ( jac, dim );

    inversa ( jac, inv_jac, dim, det_jac );

    mat_B ( B, dim, tipo_elem, nodos_elem, inv_jac, p, n, c );

    mult_mat_vector ( e, B, u, rB, cB, cB );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
