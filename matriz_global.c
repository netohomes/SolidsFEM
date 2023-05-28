//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                    /** MATRIZ DE RIGIDEZ GLOBAL **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "principal.h"
#include "memoria.h"
#include "matriz_global.h"
#include "matriz_elemental.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void ady ( struct elemental *elem, int n_elem, int nodos, int nodos_elem )
{
    /** Descripción **/
    // Obtiene la tabla de adyacencias a partir de la conectividad
    // proporcionada por GID. La conectividad es guardada en el vector
    // conec de la estructura elemental.

    /** Variables de entrada **/
    // n_elem               Num. de elementos de la malla
    // nodos                Num. de nodos de la malla
    // nodos_elem           Num. de nodos del elemento que se empleó
    //                      para mallar
    // elem[].conec         Contiene los nodos de los elementos

    /** Variables de salida **/
    // tabla_ady[].cols     Guarda el total de columnas de la matriz rala
    //                       comprimida
    // tabla_ady[].ind[]    Contiene los í­ndices de los valores no cero
    //                      de la matriz
    //                      de rigidez sin comprimir

    int i, j, k, h, a;
    int pos, no;

    tabla_ady = mem_vector_ady ( nodos, "de tabla de adyacencias");
    for ( i = 0; i < nodos; i++ ){
        tabla_ady[ i ].ind = mem_vector_int (1,"de indices de tabla_ady");
        tabla_ady[ i ].cols = 1;
    }

    for ( i = 0; i < n_elem; i++){
        for ( j = 0; j < nodos_elem; j++ ){
            for ( k = 0; k < nodos_elem; k++){
                if ( k == j ) continue;


                // Compara con todos los elementos existentes en la
                // tabla de adyacencias el valor que se quiere incor-
                // porar para no repetirlo
                pos = 1;
                no = 0;
                //a = elem[ i ].conec[ j ] - 1;
                a = elem[ i ].conec[ j ];
                for ( h = 1; h < tabla_ady[ a ].cols; h++ ){
                    if ( elem[ i ].conec[ k ] != tabla_ady[ a ].ind[ h ] ){
                        pos++;
                    }
                    else{
                        no = 1;
                    }
                }


                if ( no == 1 ) continue;

                pos++;

                // Si el valor no existe, pide memoria para el nuevo
                // elemento
                redim_vector_ady ( tabla_ady, a, pos, "adyacencias" );

                // Incorpora el nuevo elemento
                tabla_ady[ a ].ind[ pos - 1 ] = elem[i].conec[k];

                tabla_ady[ a ].cols = pos;

            }
        }
    }

    // Inicializa la primer columna de la tabla_ady que correponde
    // al núm. de nodo
    for ( i = 0; i < nodos; i++){
        tabla_ady[ i ].ind[ 0 ] = i;
    }

    /*printf ( "\n\nTABLA ADYACENCIAS \n\n" );
    for ( i = 0; i < nodos; i++){
        for ( j = 0; j < tabla_ady[ i ].cols; j++ ){
            printf ( "%d\t", tabla_ady[ i ].ind[ j ] );
            if ( j == tabla_ady[ i ].cols - 1 ) printf( "\n");
        }
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void rigid ( struct adyacencias *tabla_ady, struct elemental *elem,
             int *f_cpo_ind, double **f_cpo_val, int tipo_elem,
             int nodos_elem, int pg_elem, int dim )
{
    /** Descripción **/
    // Obtiene las matrices elementales y las ensambla en el formato
    // de renglones comprimidos

    /** Variables de entrada **/
    // tipo_elem            Tipo de elemento que se empleó en la malla
    // nodos_elem           Núm. de nodos del elemento
    // pg_elem              Núm. de puntos de Gauss para integrar de
    //                      acuerdo al tipo de elemento
    // dim                  Dimensión
    // tabla_ady[].cols     Guarda el total de columnas de la matriz
    //                      rala comprimida
    // tabla_ady[].ind[]    Contiene los índices de los valores no cero
    //                      de la matriz de rigidez sin comprimir

    /** Variables de salida **/
    // tabla_ady[].val[]    Matriz de Rigidez Global en formato
    //                      comprimido

    /** Variables locales **/
    // K_elem[][]           Matriz de rigidez elemental
    // coord_aux            Matriz con las coordenadas globales de los
    //                      nodos de un elemento
    // prop_aux[]           Vector con las propiedades de un elemento
    // f_elem[]             Vector de fuerzas debidas a las fuerzas de
    //                      cuerpo

    int i, j, m, k, r, c, pr, pc, qr, qc;
    int a, b;
    double **K_elem, **coord_aux, **D, *f_elem, **pg, **wpg, **jac, **B;
    double **BT, **Baux, **K_elemaux, **derin, **NT, *f_elem_aux, *g, E, v;

    // Memoria para las variables locales
    a         = nodos_elem;
    K_elem    = mem_matriz ( dim * a, dim * a, "K_elem." );
    D         = mem_matriz ( 6, 6, "Constitutiva." );
    f_elem    = mem_vector ( dim * a, "fuerzas de cuerpo." );
    pg        = mem_matriz ( pg_elem, dim, "de puntos de Gauss." );
    wpg       = mem_matriz ( pg_elem, dim, "de pesos de P.G." );
    jac       = mem_matriz ( dim, dim, "del Jacobiano." );
    B         = mem_matriz ( 6, a * dim, "B." );
    BT        = mem_matriz ( a * dim, 6, "BT." );
    derin     = mem_matriz ( dim, a,  "de derivadas naturales." );
    Baux      = mem_matriz ( 6, dim * a,  "Baux." );
    coord_aux = mem_matriz ( a, dim, "de coord_aux para elems.");
    K_elemaux = mem_matriz ( dim * a, dim * a, "K_elemaux." );
    det       = mem_vector ( n_elem, "determiantes de elems." );
    NT        = mem_matriz ( dim * a, dim, "NT." );
    g         = mem_vector ( dim, "peso propio." );
    f_elem_aux= mem_vector ( dim * a, "k_elem_aux." );

    //FILE *fp;
    //fp = fopen( "matriz_elemental_completa.txt", "w" );

    // Memoria para los valores de la matriz de rigidez
    // Inicializa matriz de rigidez y vector de flujos globales
    K = ( double ** ) calloc ( dim * nodos, sizeof ( double * ) );
    F = mem_vector ( dim * nodos, "fuerzas F." );
    U = mem_vector ( dim * nodos, "solucion" );
    for ( i = 0; i < nodos; i++ ){
        a  = dim * tabla_ady[ i ].cols;
        pr = i * dim;
        for ( j = 0; j < dim; j++ ){
            K[ pr + j ] = mem_vector ( a, "matriz de rigidez K." );
        }
        //for ( j = 0; j < a; j++ )   K[ i ][ j ] = 0.0;
        for ( j = 0; j < dim; j++ ) F[ i + j ]  = 0.0;
    }

    // Inicializa P.G. y pesos
    P_G ( pg, wpg, tipo_elem, dim );

    for ( m = 0; m < n_elem; m++ ){

        // Inicializa la matriz de rigidez elemental.
        inicializa_matriz ( K_elem, nodos_elem, nodos_elem );

        // Inicializa vector de fuerzas de cuerpo
        inicializa_vector ( f_elem, nodos_elem );

        // Obtiene las coordenadas de  los nodos del elemento
        for ( i = 0; i < nodos_elem; i++){
            for ( j = 0; j < dim; j++){
                a = elem[ m ].conec[ i ];
                coord_aux[ i ][ j ] = coord[ a ][ j ];
            }
        }

        // Obtiene el vector de fuerzas del elemento
        for ( i = 0; i < dim; i++ ) g[ i ] = f_cpo_val[ m ][ i ];

        /*printf( "\n\n ELEMENTO %d\n", m );
        printf( "\n Coordenadas\n");
        for ( i = 0; i < nodos_elem; i++){
            for ( j = 0; j < dimension; j++){
                printf ( "%lf\t", coord_aux[i][j]);
                if ( j == dimension - 1) printf( "\n");
            }
        }*/

        // Evalua la matriz de rigidez elemental para el elemento m
        E = elemento[ m ].E;
        v = elemento[ m ].v;

        rigid_elem ( K_elem, f_elem, D, coord_aux, B, BT, det, pg, wpg,
                     jac, derin, Baux, K_elemaux, dim, tipo_elem,
                     nodos_elem, pg_elem, NT, g, f_elem_aux, m, E, v );

        /*printf ( "\nMatriz k_elem\n " );
        for ( a = 0; a < nodos_elem * dim; a++ ){
            for ( b = 0; b < nodos_elem * dim; b++ ){
                printf ( "%8.3lf\t", K_elem[a][b] );
            if ( b == dim * nodos_elem - 1 ) printf ( "\n" );
            }
        }*/

        /*fprintf ( fp, "\nMatriz k_elem\n " );
        for ( a = 0; a < nodos_elem * dim; a++ ){
            for ( b = 0; b < nodos_elem * dim; b++ ){
                fprintf ( fp, "%8.3lf\t", K_elem[a][b] );
            if ( b == dim * nodos_elem - 1 ) fprintf ( fp, "\n" );
            }
        }*/



        /*printf ( "\n\nf_elem " );
        for ( a = 0; a < dim * nodos_elem; a++ ){
            printf ( "\n%lf", f_elem[ a ] );

        }*/

        // Ensambla las matrices elementales y vector de flujos
        for ( i = 0; i < nodos_elem; i++ ){
            for ( j = 0; j < nodos_elem; j++ ){
                a = elem[ m ].conec[ i ]; // renglón
                b = elem[ m ].conec[ j ]; // columna
                pr = a * dim;
                qr = i * dim;
                qc = j * dim;
                for ( k = 0; k < tabla_ady[ a ].cols; k++){
                    if ( tabla_ady[ a ].ind[ k ] == b ){
                        pc = k * dim;
                        for ( r = 0; r < dim; r++ ){
                            for ( c = 0; c < dim; c++ ){
                                K[ pr + r ][ pc + c ] += K_elem[ qr + r ][ qc + c ];
                            }
                        }
                        //printf ( "\nK [ %d ][ %d ] = %lf", i, j, K_elem[ i ][ j ]);

                        //printf ( "\ntabla_ady[ %d ].val[ %d ] = %lf\n", a, b, tabla_ady[ a ].val[ b ] );
                    }
                }
            }
            for ( k = 0; k < dim; k++ ){
                F[ pr + k ] += f_elem[ qr + k ];
            }
            //tabla_ady[ a ].f += f_elem[ i ];
        }
    }

    /*printf ( "\n\nMATRIZ DE RIGIDEZ COMPRIMIDA \n\n" );
    for ( i = 0; i < nodos; i++){
        pr = i * dim;
        for ( j = 0; j < dim; j++ ){
            for ( k = 0; k < dim * tabla_ady[ i ].cols; k++ ){
                printf ( "%lf\t", K[ pr + j ][ k ] );
                if ( k == dim * tabla_ady[ i ].cols -1 ) printf( "\n");
            }
        }
    }

    printf ( "\n\nVECTOR DE FUERZAS DE CUERPO \n" );
    for ( i = 0; i < nodos * dim; i++ ){
        printf ( "\n%lf", F[ i ] );
    }*/

    /*libera_matriz ( K_elem, nodos_elem );
    libera_matriz ( coord_aux, nodos_elem );
    libera_vector ( f_elem );
    libera_vector ( prop_aux );

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        libera_matriz ( pg, pg_elem );
        libera_matriz ( wpg, pg_elem );
        libera_matriz ( B, dimension );
        libera_matriz ( BT, nodos_elem );
        libera_matriz ( Baux, dimension );
        libera_matriz ( K_elemaux, nodos_elem );
        libera_matriz ( derin, dimension );
        libera_matriz ( difusion, dimension );
        libera_vector ( vect_aux );
        libera_vector ( vect_aux2 );
        libera_matriz ( jacobiano, dimension );
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void P_G ( double **pg, double **wpg, int tipo_elem, int dim )
{

	/** Descripción **/
	// Devuelve los puntos de Gauss y los pesos que deberán ser evaluados según
	// el tipo de elemento que se emplee

	/** Variables de entrada **/
	// pg[][]		Apuntador a la matriz pg
	// wpg[][]		Apuntador a la matriz wpg
	// tipo_elem	Tipo de elemento que depende del número de nodos y el número de
	//			puntos de Gauss que se emplearán

	/** Variables de salida **/
	// pg[][]		Contine los puntos de Gauss
	// wpg[][]		Contiene los pesos

    if ( dim == 3 ){

        switch ( tipo_elem ){
            case 3:
            {

                pg[0][0] = 0.25;
                pg[0][1] = 0.25;
                pg[0][2] = 0.25;
                wpg[0][0] = 1.0 / 6.0;
                wpg[0][1] = 1.0;
                wpg[0][2] = 1.0;

                break;
            }
            case 4:
            {

                pg[0][0] = - sqrt ( 3.0 ) / 3.0;
                pg[0][1] = - sqrt ( 3.0 ) / 3.0;
                pg[0][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[0][0] = 1.0;
                wpg[0][1] = 1.0;
                wpg[0][2] = 1.0;

                pg[1][0] =   sqrt ( 3.0 ) / 3.0;
                pg[1][1] = - sqrt ( 3.0 ) / 3.0;
                pg[1][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[1][0] = 1.0;
                wpg[1][1] = 1.0;
                wpg[1][2] = 1.0;

                pg[2][0] =   sqrt ( 3.0 ) / 3.0;
                pg[2][1] =   sqrt ( 3.0 ) / 3.0;
                pg[2][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[2][0] = 1.0;
                wpg[2][1] = 1.0;
                wpg[2][2] = 1.0;

                pg[3][0] = - sqrt ( 3.0 ) / 3.0;
                pg[3][1] =   sqrt ( 3.0 ) / 3.0;
                pg[3][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[3][0] = 1.0;
                wpg[3][1] = 1.0;
                wpg[3][2] = 1.0;

                pg[4][0] = - sqrt ( 3.0 ) / 3.0;
                pg[4][1] = - sqrt ( 3.0 ) / 3.0;
                pg[4][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[4][0] = 1.0;
                wpg[4][1] = 1.0;
                wpg[4][2] = 1.0;

                pg[5][0] =   sqrt ( 3.0 ) / 3.0;
                pg[5][1] = - sqrt ( 3.0 ) / 3.0;
                pg[5][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[5][0] = 1.0;
                wpg[5][1] = 1.0;
                wpg[5][2] = 1.0;

                pg[6][0] =   sqrt ( 3.0 ) / 3.0;
                pg[6][1] =   sqrt ( 3.0 ) / 3.0;
                pg[6][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[6][0] = 1.0;
                wpg[6][1] = 1.0;
                wpg[6][2] = 1.0;

                pg[7][0] = - sqrt ( 3.0 ) / 3.0;
                pg[7][1] =   sqrt ( 3.0 ) / 3.0;
                pg[7][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[7][0] = 1.0;
                wpg[7][1] = 1.0;
                wpg[7][2] = 1.0;

                break;
            }
            default:
            {
                puts ( "Error: no se pudieron determinar los puntos de Gauss. Caso no incluido. " );
                exit ( 1 );
            }
        }
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

