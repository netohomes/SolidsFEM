//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                                /** MATRICES ELEMENTALES **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "matriz_elemental.h"
#include "oper_matrices.h"
#include "memoria.h"
#include <stdio.h>
#include <stdlib.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void rigid_elem ( double **K_elem, double *f_elem, double **D, double **coord, double **B,
                  double **BT, double *det, double **pg, double **wpg, double **jac,
                  double **derin, double **Baux, double **K_elemaux, int dim, int tipo_elem,
                  int nodos_elem, int pg_elem, double **NT, double *g, double *f_elem_aux,
                  int m, double E, double v )
{
    /** Descripción **/
    // Calcula la matriz de rigidez de un elemento y su vector de fuerzas debido a las
    // fuerzas gravitatiorias.

    /** Variables de entrada **/
    // [D]              Para almacenar la ecuación constitutiva.
    // [K_elem]         Para almacenar la matriz de rigidez del elemento.
    // [K_elemaux]      Para almacenar la matriz de rigidez elemental asociada a un PG.
    // {f_elem}         Para guardar el vector de fuerzas gravitatorias.
    // [coord]          Coordenadas globales de los nodos del elemento.
    // [B] y [BT]       Para guardar la matriz B y su transpuesta
    // [Baux]           Para alamacenar el procudto de [D][B]
    // [pg] y [wpg]     Coordenas y pesoso de los PG.
    // [jac]            Para guardar la matriz jacobiana.
    // [derin]          Para guardar la matriz de derivadas naturales.
    // dim              Dimensión en la que se está trabajando.
    // tipo_elem        Tipo de elemento empleado para mallar.
    // nodos_elem       Nodos de este elemento.
    // pg_elem          Número de PG empleados para la integración Gaussiana.
    // {det}            Para almacenar los determinantes de los elementos.
    // m                Número de elemento.


    /** Variables de salida **/
    // [K_elem]         Matriz de rigidez del elemento en cuestión.
    // {f_elem}         Vector de fuerzas gravitatorias del elemento en cuestión.
    // {det}            Determinante del elemento m.

	/** Variables locales **/
	int k;
	double **inv_jac;
	double   p, n, c, wp , wn, wc;
	double   det_jac;
	int rB, cB, rD, cD, rf;
	inv_jac = mem_matriz ( dim, dim,  "del inv_jacobiano." );
    //
    int a, b;

	//FILE *fp;
	//fp = fopen( "Matriz_elem.txt", "w" );

    // Tamaño de matrices
	tam_mat ( tipo_elem, dim, nodos_elem, &rB, &cB, &rD, &cD, &rf );


    // Inicialización de variabes principales
    inicializa_vector ( f_elem, rf );
    inicializa_matriz ( K_elem, cB, cB );

	// Determina la matriz constitutiva
    constitutiva ( D, E, v, dim );

    for ( k = 0; k < pg_elem; k++ ){

        PG_k ( pg, wpg, dim, k, &p, &n, &c, &wp, &wn, &wc );

            inicializa_matriz ( inv_jac, dim, dim );
            inicializa_matriz ( B, rB, cB );
            inicializa_matriz ( BT, cB, rB );
            inicializa_matriz ( Baux, rD, cB );
            inicializa_matriz ( K_elemaux, cB, cB );
            inicializa_vector ( f_elem_aux, rf );

        mat_jacobiana ( jac, coord, p, n, c, nodos_elem, tipo_elem, dim, derin );

        det_jac = determinante ( jac, dim );

        det[ m ] = det_jac;

            //printf ( "\ndet_jacobiano = %lf\n\n", det_jac );

        inversa ( jac, inv_jac, dim, det_jac );
            /*
            printf ( "\njacobiano\n");
            for ( a = 0; a < dim; a++ ){
                for ( b = 0; b < dim; b++ ){
                    printf ( "%lf\t", jac[a][b] );
                    if ( b == dim - 1 ) printf ( "\n" );
                }
            }
            printf ( "\n\n" );

            printf ( "\ninversa_jacobiano\n");
            for ( a = 0; a < dim; a++ ){
                for ( b = 0; b < dim; b++ ){
                    printf ( "%lf\t", inv_jac[a][b] );
                    if ( b == dim - 1 ) printf ( "\n" );
                }
            }
            printf ( "\n\n" );
            */
        mat_B ( B, dim, tipo_elem, nodos_elem, inv_jac, p, n, c );
            /*fprintf ( fp, "\nPunto de Gauss: %lf\t%lf\t%lf\n", p, n, c );
            printf ( "\nB\n " );
            for ( a = 0; a <  rB; a++ ){
                for ( b = 0; b < cB; b++ ){
                    fprintf ( fp, "%lf\t", B[a][b] );
                //if ( b == cB - 1 ) fprintf ( fp, "\n" );
                }
            }*/

        transponer_mat ( B, BT, rB, cB );

            /*printf ( "\nBT\n " );
            for ( a = 0; a <  cB; a++ ){
                for ( b = 0; b < rB; b++ ){
                    printf ( "%lf\t", BT[a][b] );
                if ( b == dim - 1 ) printf ( "\n" );
                }
            }*/

        multiplica ( Baux, D, B, rD, cD, cB );

            //printf ( "\nBaux\n " );
            /*fprintf ( fp, "\n%lf\t%lf\t%lf\n", p, n, c );
            for ( a = 0; a <  rD; a++ ){
                for ( b = 0; b < cB; b++ ){
                    fprintf ( fp, "%.2lf\t",Baux[a][b] );
                //if ( b == nodos_elem - 1 ) printf ( "\n" );
                }
            }*/
            /*fprintf ( fp, "\n%lf\t%lf\t%lf\n", p, n, c );
            for ( a = 0; a <  rD; a++ ){
                for ( b = 0; b < cD; b++ ){
                    fprintf ( fp, "%.2lf\t",D[a][b] );
                //if ( b == nodos_elem - 1 ) printf ( "\n" );
                }
            }*/

        multiplica       ( K_elemaux, BT, Baux, cB, rB, cB );

        mult_mat_escalar ( K_elemaux, wc, nodos_elem, nodos_elem );

            mult_mat_escalar ( K_elemaux, det_jac, cB, cB );
            mult_mat_escalar ( K_elemaux, wp, cB, cB );
            mult_mat_escalar ( K_elemaux, wn, cB, cB );
            mult_mat_escalar ( K_elemaux, wc, cB, cB );

            /*printf ( "\nMatriz k_elemaux\n " );
            for ( a = 0; a <  dim * nodos_elem; a++ ){
                for ( b = 0; b < dim * nodos_elem; b++ ){
                    printf ( "%lf\t", K_elemaux[a][b] );
                if ( b == dim * nodos_elem - 1 ) printf ( "\n" );
                }
            }*/
            /*fprintf ( fp, "\n\n%lf\t%lf\t%lf\n", p, n, c );
            for ( a = 0; a <  dim * nodos_elem; a++ ){
                for ( b = 0; b < dim * nodos_elem; b++ ){
                    fprintf ( fp, "%lf\t", K_elemaux[a][b] );
                //if ( b == dim * nodos_elem - 1 ) printf ( "\n" );
                }
            }*/

        // Suma las matrices elementales corresponcientes a cada P.G.
        suma_mat ( K_elem, K_elemaux, cB, cB );

        fuer_elem ( f_elem_aux, NT, g, p, n, c, dim, tipo_elem, nodos_elem );

            mult_vec_escalar ( f_elem_aux, det_jac, rf );
            mult_vec_escalar ( f_elem_aux, wp, rf );
            mult_vec_escalar ( f_elem_aux, wn, rf );
            mult_vec_escalar ( f_elem_aux, wc, rf );

        suma_vec ( f_elem, f_elem_aux, rf );



        // Imprime en pantalla f_elem para cada P.G.
        /*printf ( "\n\nf_elem " );
        for ( a = 0; a < dim * nodos_elem; a++ ){
            printf( "\n%lf", f_elem[ a ] );
        }*/

    }

    //fclose( fp );

    libera_matriz ( inv_jac, dim );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void fuer_elem ( double *f_elem_aux, double **NT, double *b, double p, double n, double c, int dim,
                 int tipo_elem, int nodos_elem )
{
    /** Descripción **/
    // Calcula el vecto de fuerzas debido a la gravedad

    int rNT, cNT;

    rNT = dim * nodos_elem;
    cNT = dim;

    mat_N ( NT, p, n, c, dim, tipo_elem, nodos_elem, rNT, cNT );

    mult_mat_vector ( f_elem_aux, NT, b, rNT, cNT, dim );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void mat_N ( double **NT, double p, double n, double c, int dim, int tipo_elem, int nodos_elem,
            int rNT, int cNT )
{
    /** Descripción **/
    // Calcula [N] para un PG dado.

    int i, j;
    int a;
    double Ni;

    inicializa_matriz ( NT, rNT, cNT );

    for ( i = 0; i < nodos_elem; i++ ){

        Ni = N ( tipo_elem, i, dim, p, n, c );

        for ( j = 0; j < dim; j++ ){

            a = dim * i;
            NT[ a + j ][ j ] = Ni;
            NT[ a + j ][ j ] = Ni;

        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void mat_B ( double **B, int dim, int tipo_elem, int nodos_elem, double **inv_jac, double p,
             double n, double c )
{
    /** Descripción **/
    // Calcula la matriz B para un dado P.G.

    int j;
    double *v1, *v2;
    v1 = mem_vector ( nodos_elem, "vector aux. v1."   );
	v2 = mem_vector ( nodos_elem, "vector aux. v2."   );

    for ( j = 0; j < nodos_elem; j++ ){

        inicializa_vector ( v1, dim );
        inicializa_vector ( v2, dim );

        deriv_Ni ( v1, p, n, c, tipo_elem, dim, j );

            /*printf ( "\n\nvect_aux" );
            for ( a = 0; a <  dim; a++ ){
                printf ( "\n%lf", v1[a] );
            }*/

        mult_mat_vector ( v2, inv_jac, v1, dim, dim, dim );

            /*printf ( "\n\nvect_aux2" );
            for ( a = 0; a <  dim; a++ ){
                printf ( "\n%lf", vect_aux2[a] );
            }*/

        sum_B ( B, v2, dim, j );

    }

    /*printf ( "\nMatriz B\n " );
    for ( a = 0; a <  dim; a++ ){
        for ( b = 0; b < nodos_elem; b++ ){
            printf ( "%lf\t", B[a][b] );
            if ( b == nodos_elem - 1 ) printf ( "\n" );
        }
    }*/

    libera_vector ( v1 );
    libera_vector ( v2 );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void tam_mat ( int tipo_elem, int dim, int nodos_elem, int *rB, int *cB, int *rD, int *cD, int *rf )
{
    /** Descripción **/
    // Determina los tamaños de las matrices empleadas para calcular la matriz de rigiddez elemental.
    // El sufijo r hace referencia a los renglones y el sufijo c a las columnas.

    if ( dim == 3 ){

        // Para [B]
        *rB = 6;
        *cB = nodos_elem * dim;

        // Para [D]
        *rD = 6;
        *cD = 6;

        // Para {f_elem}
        *rf = *cB;

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void sum_B ( double **B, double *v2, int dim, int j )
{
    /** Descripción **/
    // Incorpora loos valores de la matriz B.
    // B debió haber sido inicializadas antes.

    int p;

    if ( dim == 3 ){

        p = j * dim;
        B[ 0 ][ p + 0 ] += v2[ 0 ];
        B[ 1 ][ p + 1 ] += v2[ 1 ];
        B[ 2 ][ p + 2 ] += v2[ 2 ];
        B[ 3 ][ p + 0 ] += v2[ 1 ];
        B[ 3 ][ p + 1 ] += v2[ 0 ];
        B[ 4 ][ p + 0 ] += v2[ 2 ];
        B[ 4 ][ p + 2 ] += v2[ 0 ];
        B[ 5 ][ p + 1 ] += v2[ 2 ];
        B[ 5 ][ p + 2 ] += v2[ 1 ];

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void mat_jacobiana ( double **jacobiano, double **coord, double p, double n, double c, int nodos_elem,
                     int tipo_elem, int dim, double **derin )
{
	/** Descripción **/
	// Evalúa el jacobiano con las coordenadas globales de los nodos
	// y un punto de Gauss que se empleará para la integración.

	/** Variables de entrada **/
	// dim	Dimensión del problema
	// nodos_elem	Número de nodos del elemento
	// coord[][]	Coordenadas globales de los nodos del elemento
	// p, n, c  	Coordenada del P.G.
	// tipo_elem	Tipo de elemento empleado

	/** Variables de salida **/
	// jacobiano	Jacobiano evaluado en P.G. ( p , n )

	/** Variables locales **/
	// derin		Derivadas naturales evaluadas en un punto de Gauss
	//			    en P.G. ( p , n )

	//int a, b;
	//double **derin;


	//derin = mem_matriz ( dim, nodos_elem,  "de derivadas naturales" );

	inicializa_matriz ( jacobiano, dim, dim );

	deriv_nat ( derin, p, n, c, tipo_elem );
    //printf( "Despues de fderin \n");
	/*printf ( "\nderin\n");
	for ( a = 0; a < dim; a++ ){
		for ( b = 0; b < nodos_elem; b++ ){
			printf ( "%lf\t", derin[a][b] );
			if ( b == nodos_elem - 1 ) printf ( "\n" );
		}
	}
	printf ( "\n\n" );*/

	multiplica ( jacobiano, derin, coord, dim, nodos_elem, dim );

	/*printf ( "\ncoord\n");
	for ( a = 0; a < nodos_elem; a++ ){
		for ( b = 0; b < dim; b++ ){
			printf ( "%lf\t", coord[a][b] );
			if ( b == dim - 1 ) printf ( "\n" );
		}
	}
	printf ( "\n\n" );*/


}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void deriv_nat ( double **derin, double p, double n, double c, int tipo_elem )
{
	/** Descripción **/
	// Evalúa las derivadas naturales en los puntos de Gauss de acuerdo al tipo de
	// elemento que se esté empleando

	/** Variables de entrada */
	// p, n, c 	Coordenada del P.G.

	/**Variables de salida **/
	// derin	Derivadas naturales evaluadas en el P.G. ( p, n, c )

	/*  c
		|    n
		|   /
		|  /
		| /
		|/
		|-------------- p
	*/

    switch ( tipo_elem ){
        case 3:
        {
            derin[0][0] = -1.0;
            derin[0][1] =  1.0;
            derin[0][2] =  0.0;
            derin[0][3] =  0.0;

            derin[1][0] = -1.0;
            derin[1][1] =  0.0;
            derin[1][2] =  1.0;
            derin[1][3] =  0.0;

            derin[2][0] = -1.0;
            derin[2][1] =  0.0;
            derin[2][2] =  0.0;
            derin[2][3] =  1.0;

            break;
        }
        case 4:
        {
            derin[0][0] = -1.0/8.0 * ( 1 - n ) * ( 1 - c );
            derin[0][1] =  1.0/8.0 * ( 1 - n ) * ( 1 - c );
            derin[0][2] =  1.0/8.0 * ( 1 + n ) * ( 1 - c );
            derin[0][3] = -1.0/8.0 * ( 1 + n ) * ( 1 - c );
            derin[0][4] = -1.0/8.0 * ( 1 - n ) * ( 1 + c );
            derin[0][5] =  1.0/8.0 * ( 1 - n ) * ( 1 + c );
            derin[0][6] =  1.0/8.0 * ( 1 + n ) * ( 1 + c );
            derin[0][7] = -1.0/8.0 * ( 1 + n ) * ( 1 + c );

            derin[1][0] = -1.0/8.0 * ( 1 - p ) * ( 1 - c );
            derin[1][1] = -1.0/8.0 * ( 1 + p ) * ( 1 - c );
            derin[1][2] =  1.0/8.0 * ( 1 + p ) * ( 1 - c );
            derin[1][3] =  1.0/8.0 * ( 1 - p ) * ( 1 - c );
            derin[1][4] = -1.0/8.0 * ( 1 - p ) * ( 1 + c );
            derin[1][5] = -1.0/8.0 * ( 1 + p ) * ( 1 + c );
            derin[1][6] =  1.0/8.0 * ( 1 + p ) * ( 1 + c );
            derin[1][7] =  1.0/8.0 * ( 1 - p ) * ( 1 + c );

            derin[2][0] = -1.0/8.0 * ( 1 - p ) * ( 1 - n );
            derin[2][1] = -1.0/8.0 * ( 1 + p ) * ( 1 - n );
            derin[2][2] = -1.0/8.0 * ( 1 + p ) * ( 1 + n );
            derin[2][3] = -1.0/8.0 * ( 1 - p ) * ( 1 + n );
            derin[2][4] =  1.0/8.0 * ( 1 - p ) * ( 1 - n );
            derin[2][5] =  1.0/8.0 * ( 1 + p ) * ( 1 - n );
            derin[2][6] =  1.0/8.0 * ( 1 + p ) * ( 1 + n );
            derin[2][7] =  1.0/8.0 * ( 1 - p ) * ( 1 + n );

            break;
        }
        default:
        {
            printf ( "\nError: no se pudieron calcular las derivadas naturales. " );
            printf ( "Caso no incluido." );
            exit   ( 1 );
        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
double determinante ( double **A, int dim )
{
	/** Descripción **/
	// Calcula el determinante de una matriz [A] de tamaño 3 x 3.

	double det;

	switch ( dim ){
		case 3:
		{
		    det =   A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*( A[1][0]*A[2][2]
                  - A[1][2]*A[2][0] ) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
            break;
		}
		default:
		{
			puts ( "Error: no se pudo calcular el determinante. Caso no incluido." );
			exit ( 1 );
		}
	}

    if ( det == 0 ){
        puts ( "Determinante del jacobiano = 0" );
        exit ( 1 );
    }

	return det;

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void inversa ( double **A, double **inversa_A, int n, double det )
{
	/** Descripción **/
	// Calcula la inversa de matrices de 3 x 3

	/** Variables de entrada **/
	// A			Matriz a invertir
	// n			Tamaño de la matriz [A]

	/** Variables de salida **/
	// inversa_A	Inversa de [A]

	switch ( n ) {
		case 3:
		{
            inversa_A[0][0] =  1.0/det * ( A[1][1]*A[2][2] - A[1][2]*A[2][1] );
            inversa_A[0][1] = -1.0/det * ( A[0][1]*A[2][2] - A[0][2]*A[2][1] );
            inversa_A[0][2] =  1.0/det * ( A[0][1]*A[1][2] - A[0][2]*A[1][1] );

            inversa_A[1][0] = -1.0/det * ( A[1][0]*A[2][2] - A[1][2]*A[2][0] );
            inversa_A[1][1] =  1.0/det * ( A[0][0]*A[2][2] - A[0][2]*A[2][0] );
            inversa_A[1][2] = -1.0/det * ( A[0][0]*A[1][2] - A[0][2]*A[1][0] );

            inversa_A[2][0] =  1.0/det * ( A[1][0]*A[2][1] - A[1][1]*A[2][0] );
            inversa_A[2][1] = -1.0/det * ( A[0][0]*A[2][1] - A[0][1]*A[2][0] );
            inversa_A[2][2] =  1.0/det * ( A[0][0]*A[1][1] - A[0][1]*A[1][0] );

		    break;
		}
		default:
		{
			puts ( "Error: no se pudo calcular la inversa. Caso no incluido." );
			exit ( 1 );
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void deriv_Ni( double *vect_aux, double p, double n, double c, int tipo_elem, int dim, int j )
{
    /** Descripción **/
    // Calcula las derivadas de la función de forma Ni respecto a las coordenadas del espacio
    // normalizado (p,n,c)

    switch ( tipo_elem ) {

        case 3:
        {
            if ( j == 0 ){
                vect_aux[0] = -1.0;
                vect_aux[1] = -1.0;
                vect_aux[2] = -1.0;
            }
            else if ( j == 1 ){
                vect_aux[0] = 1.0;
                vect_aux[1] = 0.0;
                vect_aux[2] = 0.0;
            }
            else if ( j == 2){
                vect_aux[0] = 0.0;
                vect_aux[1] = 1.0;
                vect_aux[2] = 0.0;
            }
            else{
                vect_aux[0] = 0.0;
                vect_aux[1] = 0.0;
                vect_aux[2] = 1.0;
            }
            break;
        }

        case 4:
        {
            if ( j == 0){
                vect_aux[0] = -1.0/8.0 * ( 1.0 - n ) * ( 1.0 - c );
                vect_aux[1] = -1.0/8.0 * ( 1.0 - p ) * ( 1.0 - c );
                vect_aux[2] = -1.0/8.0 * ( 1.0 - p ) * ( 1.0 - n );
            }
            else if ( j == 1 ){
                vect_aux[0] =  1.0/8.0 * ( 1.0 - n ) * ( 1.0 - c );
                vect_aux[1] = -1.0/8.0 * ( 1.0 + p ) * ( 1.0 - c );
                vect_aux[2] = -1.0/8.0 * ( 1.0 + p ) * ( 1.0 - n );
            }
            else if ( j == 2 ){
                vect_aux[0] =  1.0/8.0 * ( 1.0 + n ) * ( 1.0 - c );
                vect_aux[1] =  1.0/8.0 * ( 1.0 + p ) * ( 1.0 - c );
                vect_aux[2] = -1.0/8.0 * ( 1.0 + p ) * ( 1.0 + n );
            }
            else if ( j == 3 ){
                vect_aux[0] = -1.0/8.0 * ( 1.0 + n ) * ( 1.0 - c );
                vect_aux[1] =  1.0/8.0 * ( 1.0 - p ) * ( 1.0 - c );
                vect_aux[2] = -1.0/8.0 * ( 1.0 - p ) * ( 1.0 + n );
            }
            else if ( j == 4 ){
                vect_aux[0] = -1.0/8.0 * ( 1.0 - n ) * ( 1.0 + c );
                vect_aux[1] = -1.0/8.0 * ( 1.0 - p ) * ( 1.0 + c );
                vect_aux[2] =  1.0/8.0 * ( 1.0 - p ) * ( 1.0 - n );
            }
            else if ( j == 5 ){
                vect_aux[0] =  1.0/8.0 * ( 1.0 - n ) * ( 1.0 + c );
                vect_aux[1] = -1.0/8.0 * ( 1.0 + p ) * ( 1.0 + c );
                vect_aux[2] =  1.0/8.0 * ( 1.0 + p ) * ( 1.0 - n );
            }
            else if ( j == 6 ){
                vect_aux[0] =  1.0/8.0 * ( 1.0 + n ) * ( 1.0 + c );
                vect_aux[1] =  1.0/8.0 * ( 1.0 + p ) * ( 1.0 + c );
                vect_aux[2] =  1.0/8.0 * ( 1.0 + p ) * ( 1.0 + n );
            }
            else if ( j == 7 ){
                vect_aux[0] = -1.0/8.0 * ( 1.0 + n ) * ( 1.0 + c );
                vect_aux[1] =  1.0/8.0 * ( 1.0 - p ) * ( 1.0 + c );
                vect_aux[2] =  1.0/8.0 * ( 1.0 - p ) * ( 1.0 + n );
            }

            break;
        }

        default:
        {
            printf ( "\nError: no se pudo calcular el vec_aux para B. Caso no incluido" );
            exit ( 1 );
        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void PG_k ( double **pg, double **wpg, int dim, int k, double *p, double *n, double *c,
            double *wp,  double *wn,   double *wc )
{
    /** Descripción **/
    // Establece las coordenadas y pesos correspondientes de un P.G.

    if ( dim == 3 ){

            *p  = pg [ k ][ 0 ];
            *n  = pg [ k ][ 1 ];
            *c  = pg [ k ][ 2 ];
            *wp = wpg[ k ][ 0 ];
            *wn = wpg[ k ][ 1 ];
            *wc = wpg[ k ][ 2 ];

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void constitutiva ( double **A, double E, double v, int dim )
{
    /** Descripción **/
    // Calcula la ecuación constitutiva de acuerdo al tipo de problema

    double h, u;
    //int i, j;

    if ( dim == 3 ){

            h = E * v / ( ( 1.0 + v ) * ( 1.0 - 2.0 * v ) );
            u = E / ( 2.0 * ( 1 + v ) );

            A[ 0 ][ 0 ] = h + 2.0 * u;
            A[ 0 ][ 1 ] = h;
            A[ 0 ][ 2 ] = h;
            A[ 0 ][ 3 ] = 0.0;
            A[ 0 ][ 4 ] = 0.0;
            A[ 0 ][ 5 ] = 0.0;

            A[ 1 ][ 0 ] = h;
            A[ 1 ][ 1 ] = h + 2.0 * u;
            A[ 1 ][ 2 ] = h;
            A[ 1 ][ 3 ] = 0.0;
            A[ 1 ][ 4 ] = 0.0;
            A[ 1 ][ 5 ] = 0.0;

            A[ 2 ][ 0 ] = h;
            A[ 2 ][ 1 ] = h;
            A[ 2 ][ 2 ] = h + 2.0 * u;
            A[ 2 ][ 3 ] = 0.0;
            A[ 2 ][ 4 ] = 0.0;
            A[ 2 ][ 5 ] = 0.0;

            A[ 3 ][ 0 ] = 0.0;
            A[ 3 ][ 1 ] = 0.0;
            A[ 3 ][ 2 ] = 0.0;
            A[ 3 ][ 3 ] = u;
            A[ 3 ][ 4 ] = 0.0;
            A[ 3 ][ 5 ] = 0.0;

            A[ 4 ][ 0 ] = 0.0;
            A[ 4 ][ 1 ] = 0.0;
            A[ 4 ][ 2 ] = 0.0;
            A[ 4 ][ 3 ] = 0.0;
            A[ 4 ][ 4 ] = u;
            A[ 4 ][ 5 ] = 0.0;

            A[ 5 ][ 0 ] = 0.0;
            A[ 5 ][ 1 ] = 0.0;
            A[ 5 ][ 2 ] = 0.0;
            A[ 5 ][ 3 ] = 0.0;
            A[ 5 ][ 4 ] = 0.0;
            A[ 5 ][ 5 ] = u;

    }

    /*printf ( "MATRIZ CONSTITUTIVA\n" );
    for ( i = 0; i < 6; i++ ){
        for ( j = 0; j < 6; j++ ){
            printf ( "%.3lf\t", A[ i ][ j ] );
            if ( j == 5 ) printf ( "\n" );
        }
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
double N ( int tipo_elem, int i, int dim, double p, double n, double c )
{
    /** Descripción **/
    // Evalúa las funciones de forma en un P.G.

    double N;


    if ( dim == 3 ){

        switch ( tipo_elem ){
            case 3:
            {
                if ( i == 0 ){
                    N = 1 - p - n - c;
                }
                else if ( i == 1 ){
                    N = p;
                }
                else if ( i == 2 ){
                    N = n;
                }
                else{
                    N = c;
                }
                break;
            }
            case 4:
            {
                if ( i == 0){
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 - n ) * ( 1 - c );
                }
                else if ( i == 1 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 - n ) * ( 1 - c );
                }
                else if ( i == 2 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 + n ) * ( 1 - c );
                }
                else if ( i == 3 ){
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 + n ) * ( 1 - c );
                }
                else if ( i == 4 ){
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 - n ) * ( 1 + c );
                }
                else if ( i == 5 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 - n ) * ( 1 + c );
                }
                else if ( i == 6 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 + n ) * ( 1 + c );
                }
                else{
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 + n ) * ( 1 + c );
                }
                break;
            }
            default:
            {
                puts ( "Error: no se pudo evaluar Ni en el P.G. Caso no incluido." );
                exit ( 1 );
            }
        }
    }

    return N;

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
