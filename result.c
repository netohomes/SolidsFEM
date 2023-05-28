//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                  /** ESCRIBE RESULTADOS PARA GID **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "result.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void result( char *r, int n_elem, int pg_elem, int tipo_elem,  int nodos,
             double *U, double **ESF, int dim )
{
    /** Descripción **/
    // Escribe el archivo de resultados para que sea leído por GID.

    /** Variables de entrada **/
    // r                Dirección al archivo con el que se está
    //                  trabajando, sin extensión.
    // U                Desplazmientos de los nodos.
    // pg_elem          Puntos de Gauss del elemento que se empleó
    //                  para mallar.
    // nodos            Número total de nodos de la malla.
    // tipo_elem        Tipo de elemento que se empleó para mallar.

    /** Variables locales **/
    // ruta             Dirección al archivo con el que se está trabajando,
    //                  con extensión .post.res .
    // nom              Tipo de elemento
    // a               Núm. de componentes del vector de esfuerzos.

    int i, j, b;
    char *ruta, *nom;
    FILE *ft;

    ruta = ( char *) calloc ( 256, sizeof ( char ) );

    // Tipo de elemento
    if ( tipo_elem == 3 ){
        nom = "Tetrahedra";
    }
    if ( tipo_elem == 4 ){
        nom = "Hexahedra";
    }

    strcpy ( ruta, r );
    strcat ( ruta, ".post.res" );

    ft = fopen ( ruta , "w" );

    //====================================================================

    fprintf( ft, "Gid Post Results File 1.0\n\n" );
    fprintf ( ft, "GaussPoints BoardGaussInternal ElemType %s\n", nom );
    fprintf ( ft, "Number Of Gauss Points: %d\n", pg_elem );
    fprintf ( ft, "Natural Coordinates: internal\n" );
    fprintf ( ft, "End Gausspoints\n\n" );

    //====================================================================

    fprintf ( ft, "Result Desplazamiento Menu 0 Vector OnNodes\n" );
    fprintf ( ft, "Componentnames Dx, Dy, Dz\n" );
    fprintf ( ft, "Values\n" );
    for ( i = 0; i < nodos; i++ ){
        b = i * dim;
        fprintf ( ft, "%d", i + 1 );
        for ( j = 0; j < dim; j++ ){
            fprintf ( ft, "\t%lf", U[ b + j ] );
        }
        fprintf ( ft, "\n" );
    }
    fprintf( ft, "End Values\n\n" );

    //====================================================================

    fprintf ( ft, "Result EsfPG_A Menu 0 Vector OnGaussPoints" );
    fprintf ( ft, " BoardGaussInternal\n" );
    fprintf ( ft, "Componentnames Gx, Gy, Gz\n" );
    fprintf( ft, "Values\n");
    for ( i = 0; i < n_elem; i++ ){
        b = pg_elem * i;
        for ( j = 0; j < pg_elem; j++ ){
            if ( j == 0 ){
                fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1,
                         ESF[ b + j ][ 0 ], ESF[ b + j ][ 1 ],
                         ESF[ b + j ][ 2 ] );
                         /*
                         ESF[ b + j ][ 3 ],
                         ESF[ b + j ][ 4 ], ESF[ b + j ][ 5 ] );*/
            }
            else{
                fprintf( ft, "\t%lf\t%lf\t%lf\n",
                         ESF[ b + j ][ 0 ], ESF[ b + j ][ 1 ],
                         ESF[ b + j ][ 2 ] );
                         /*ESF[ b + j ][ 3 ],
                         ESF[ b + j ][ 4 ], ESF[ b + j ][ 5 ] );*/
            }
        }
    }
    fprintf( ft, "End Values\n\n" );

    //====================================================================

    fclose ( ft );

    // Libera nombre de archivo
    free ( ruta );

}
