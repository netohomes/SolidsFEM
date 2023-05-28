//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                    /** CONDICIONES DE CONTORNO **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "principal.h"
#include "memoria.h"
#include "cond_cont.h"
#include <math.h>
#include <stdio.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

void cond_cont ( double **K, double *F, struct adyacencias *tabla_ady,
                 double **coord, int nodos, int tipo_elem, int nodos_elem,
                 int FP , int  *f_punt_ind, double **f_punt_val,
                 int FPP, int **f_pre_ind , double  *f_pre_val ,
                 int FS , int **f_sup_ind , double **f_sup_val ,
                 int DL , int  *d_lin_ind , double **d_lin_val ,
                 int DS , int **d_sup_ind , double **d_sup_val , int dim)
{
    /** Descripción **/
    // Impone las condiciones de contorno: Dirichlet y Neumann

    //int i;

    neumann ( coord, F, FP, FPP, FS, f_punt_ind, f_punt_val, f_pre_ind,
              f_pre_val, f_sup_ind, f_sup_val, tipo_elem, nodos_elem, dim );

    dirichlet ( K, F, tabla_ady, DL, DS, d_lin_ind, d_lin_val, d_sup_ind,
                d_sup_val, nodos, tipo_elem, dim );

    /*printf ( "\n\nF" );
    for ( i = 0; i < nodos * dim; i++ ){
        printf ( "\n%lf", F[ i ] );
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void dirichlet ( double **K, double *F, struct adyacencias *tabla_ady,
                 int DL, int DS, int *d_lin_ind, double **d_lin_val,
                 int **d_sup_ind, double **d_sup_val, int nodos,
                 int tipo_elem, int dim )
{
     /** Descripción **/
    // Impone las condiciones de contorno tipo dirichlet. En este caso
    // son los desplazamientos.

    if ( DL > 0 ){
        dirich_DL ( K, F, tabla_ady, d_lin_ind, d_lin_val, DL, nodos, dim );
    }

    if ( DS > 0 ){
        dirich_DS ( K, F, tabla_ady, d_sup_ind, d_sup_val, DS, nodos,
                    tipo_elem, dim );
    }


}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void dirich_DL ( double **K, double *F, struct adyacencias *tabla_ady,
                 int *d_lin_ind, double **d_lin_val, int DL, int nodos,
                 int dim )
{
    /** Descripción **/
    // Impone los desplazamientos fijados sobre líneas.

    /** Variables de entrada **/
    // {d_lin_ind}          Índides de desplazamientos en línea.
    // [d_lin_val]          Valores de desplazamientos en línea.
    // DL                   Número de nodos con esta condición.
    // [K]                  Matriz de rigidez global en formato de
    //                      renglón comprimido.
    // {F}                  Vector global de fuerzas.
    // [tabla_ady]          Tabla de adyacencias.
    // nodos                Número de nodos de la malla.
    // dim                  Dimensión.

    /** Variables locales **/
    int i, j, k, m, n, a, b, c;
    double aux;

    // Suma la contribución de las restricciones al vector de fuerzas
    for ( i = 0; i < DL; i++ ){
        a = d_lin_ind[ i ]; // núm. de nodo
        for ( j = 0; j < nodos; j++ ){
            c = j * dim;
            for ( k = 0; k < tabla_ady[ j ].cols; k++ ){
                if ( tabla_ady[ j ].ind[ k ] == a ){
                    b = k * dim;
                    for ( m = 0; m < dim; m++ ){
                        for ( n = 0; n < dim; n++ ){
                            aux = K[ c + m ][ b + n ] *
                                  d_lin_val[ i ][ n ];
                            F[ c + n ] += aux;
                            K[ c + m ][ b + n ] = 0.0;
                        }
                    }
                }
            }
        }
    }
    // Pone 0 y 1 en los renglones que se les impuso una restricción.
    // También coloca los deplazamientos impuestos en F.
    for ( i = 0; i < DL; i++ ){
        a = d_lin_ind[ i ]; // núm. de nodo
        for ( j = 0; j < nodos; j++ ){
            if ( j == a ){
                c = j * dim;
                for ( m = 0; m < dim; m++ ) F[ c + m ] = d_lin_val[ i ][ m ];
                for ( k = 0; k < tabla_ady[ j ].cols; k++ ){
                    b = k * dim;
                    if ( k == 0 ){
                        for ( m = 0; m < dim; m++ ){
                            for ( n = 0; n < dim; n++ ){
                                if ( m == n ){
                                    K[ c + m ][ b + n ] = 1.0;
                                }
                                else{
                                    K[ c + m ][ b + n ] = 0.0;
                                }
                            }
                        }
                    }
                    else{
                        for ( m = 0; m < dim; m++ ){
                            for ( n = 0; n < dim; n++ ){
                                 K[ c + m ][ b + n ] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void dirich_DS ( double **K, double *F, struct adyacencias *tabla_ady,
                 int **d_sup_ind, double **d_sup_val, int DS, int nodos,
                 int tipo_elem, int dim )
{
    /** Descripción **/
    // Impone los desplazamientos fijados sobre líneas.

    /** Variables de entrada **/
    // {d_sup_ind}          Índides de desplazamientos de superficie.
    // [d_sup_val]          Valores de desplazamientos de superficie.
    // DS                   Número de caras con esta condición.
    // [K]                  Matriz de rigidez global en formato de
    //                      renglón comprimido.
    // {F}                  Vector global de fuerzas.
    // [tabla_ady]          Tabla de adyacencias.
    // nodos                Número de nodos de la malla.
    // dim                  Dimensión.
    // tipo_elem            Tipo de elemento empleado para mallar.

    /** Variables locales **/
    int i, j, k, m, n, p, a, b, c, nc;
    double aux;

    if ( tipo_elem == 3 ) nc = 3;
    if ( tipo_elem == 4 ) nc = 4;

    // Suma la contribución de las restricciones al vector de fuerzas
    for ( i = 0; i < DS; i++ ){
        for ( p = 0; p < nc; p++ ){
            a = d_sup_ind[ i ][ p ]; // núm. de nodo
            for ( j = 0; j < nodos; j++ ){
                c = j * dim;
                for ( k = 0; k < tabla_ady[ j ].cols; k++ ){
                    if ( tabla_ady[ j ].ind[ k ] == a ){
                        b = k * dim;
                        for ( m = 0; m < dim; m++ ){
                            for ( n = 0; n < dim; n++ ){
                                aux = K[ c + m ][ b + n ] *
                                      d_sup_val[ i ][ n ];
                                F[ c + n ] += aux;
                                K[ c + m ][ b + n ] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
    // Pone 0 y 1 en los renglones que se les impuso una restricción.
    // También coloca los deplazamientos impuestos en F.
    for ( i = 0; i < DS; i++ ){
        for ( p = 0; p < nc; p++ ){
            a = d_sup_ind[ i ][ p ]; // núm. de nodo
            for ( j = 0; j < nodos; j++ ){
                if ( j == a ){
                    c = j * dim;
                    for ( m = 0; m < dim; m++ ) F[ c + m ] = d_sup_val[ i ][ m ];
                    for ( k = 0; k < tabla_ady[ j ].cols; k++ ){
                        b = k * dim;
                        if ( k == 0 ){
                            for ( m = 0; m < dim; m++ ){
                                for ( n = 0; n < dim; n++ ){
                                    if ( m == n ){
                                        K[ c + m ][ b + n ] = 1.0;
                                    }
                                    else{
                                        K[ c + m ][ b + n ] = 0.0;
                                    }
                                }
                            }
                        }
                        else{
                            for ( m = 0; m < dim; m++ ){
                                for ( n = 0; n < dim; n++ ){
                                     K[ c + m ][ b + n ] = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void neumann ( double **coord, double *F, int FP, int FPP, int FS,
               int *f_punt_ind, double **f_punt_val, int **f_pre_ind,
               double *f_pre_val, int **f_sup_ind, double **f_sup_val,
               int tipo_elem, int nodos_elem, int dim )
{
    /** Descripción **/
    // Impone las condiones tipo Neumann en el vector de fuerzas F

    /** Variables de entrada **/

    /** Variables locales **/

    int i;

    // Fuerzas puntuales
    if ( FP > 0 ){
        neu_FP ( F, f_punt_ind, f_punt_val, FP, dim );
    }

    // Fuerzas de presión
    if ( FPP > 0 ){
        neu_FPP ( F, f_pre_ind, f_pre_val, FPP, coord, tipo_elem, nodos_elem,
                  dim );
    }

    // Fuerzas de superficie
    if ( FS > 0 ){
        neu_FS ( F, f_sup_ind, f_sup_val, FS, coord, tipo_elem, nodos_elem,
                 dim );
    }

    printf ( "\n\nVECTOR DE FUERZASS CON CC NEUMANN." );

    for ( i = 0; i < dim * nodos; i++ ){
        printf ( "\n %d\t %lf", i, F[ i ] );
    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void neu_FP ( double *F, int *f_punt_ind, double **f_punt_val, int FP,
              int dim )
{
    /** Descripción **/
    // Impone las condiciones tipo neumann aplicadas en nodos o puntos

    /** Variables de entrada **/
    // {F}          Vector de fuerzas
    // FP           Número de fzas. puntuales-
    // {f_punt_ind} Puntos en los que están aplicadas las fzas. puntuales.
    // {f_punt_val} Valores de las fzas. puntuales.
    // dim          Dimensión

    int i, j, a;

    for ( i = 0; i < FP; i++ ){

        a = f_punt_ind[ i ] * dim;

        for ( j = 0; j < dim; j++ ){

            F[ a + j ] += f_punt_val[ i ][ j ];

        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void neu_FPP ( double *F, int **f_pre_ind, double *f_pre_val, int FPP,
               double **coord, int tipo_elem, int nodos_elem, int dim )
{
    /** Descripción **/
    // Impone las fuerzas de presión al vector de fuerzas. Si es positiva
    // empuja la cara y si es negativa la jala.

    /** Variables de entrada **/
    // {F}          Vector de fuerzas
    // FP           Número de fzas. puntuales-
    // {f_punt_ind} Puntos en los que están aplicadas las fzas. puntuales.
    // {f_punt_val} Valores de las fzas. puntuales.
    // tipo_elem    Tipo de elemento de la malla.
    // [coord]      Coordenadas de todos los nodos.
    // dim          Dimensión.


    /** Variables locales **/

    if ( tipo_elem == 3 ) FPP_3 ( F, f_pre_ind, f_pre_val, FPP, coord,
                                  dim );

    if ( tipo_elem == 4 ) FPP_4 ( F, f_pre_ind, f_pre_val, FPP, coord,
                                  dim );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void FPP_3 ( double *F, int **f_pre_ind, double *f_pre_val, int FPP,
             double **coord, int dim )
{
    /** Descripción **/
    // Impone las fuerzas de presión al vector de fuerzas. Si es positiva
    // empuja la cara y si es negativa la jala. Esto para elementos tipo 3
    // (tetraedros) en 3D.

    /** Variables de entrada **/
    // {F}          Vector de fuerzas
    // FP           Número de fzas. de presión.
    // {f_pre_ind}  Caras en las que está aplicada la fuerzas.
    // {f_punt_val} Valores de las fuerzas de presión.
    // [coord]      Coordenadas de todos los nodos.
    // dim          Dimensión.

    /** Variables locales **/
    int i, a, b, c, j;
    double ax, ay, az, bx, by, bz, a_b, Ma, Mb, Pab, h, area, nx, ny, nz;
    double Mn, f, fx, fy, fz;

    for ( i = 0; i < FPP; i++ ){

        a = f_pre_ind [ i ][ 0 ];
        b = f_pre_ind [ i ][ 1 ];
        c = f_pre_ind [ i ][ 2 ];

        ax = coord[ b ][ 0 ] - coord[ a ][ 0 ];
        ay = coord[ b ][ 1 ] - coord[ a ][ 1 ];
        az = coord[ b ][ 2 ] - coord[ a ][ 2 ];

        bx = coord[ c ][ 0 ] - coord[ a ][ 0 ];
        by = coord[ c ][ 1 ] - coord[ a ][ 1 ];
        bz = coord[ c ][ 2 ] - coord[ a ][ 2 ];

        //Producto punto a.b
        a_b = ax * bx + ay * by + az * bz;

        // Magnitud de a y b
        Ma = sqrt ( ax * ax + ay * ay + az * az );
        Mb = sqrt ( bx * bx + by * by + bz * bz );

        // Proyección de a sobre b
        Pab = a_b / Mb;

        // Altura
        h = sqrt ( Ma * Ma - Pab * Pab );

        // Area del triángulo
        area = 0.5 * Mb * h;

        // Fuerza total
        f = f_pre_val[ i ] * area;

        // Producto cruz a x b
        nx = ay * bz - az * by;
        ny = az * bx - ax * bz;
        nz = ax * by - ay * bx;

        // Magnitud del vector normal
        Mn = sqrt ( nx * nx + ny * ny + nz * nz );

        // Vectores unitarios en x, y, z
        nx /= Mn;
        ny /= Mn;
        nz /= Mn;

        // Componentes del vector f
        fx = f * nx;
        fy = f * ny;
        fz = f * nz;

        // Suma de componente al vector F
        for ( j = 0; j < 3; j++ ){

            a = f_pre_ind[ i ][ j ];

            F[ a * dim + 0 ] += fx / 3.0;
            F[ a * dim + 1 ] += fy / 3.0;
            F[ a * dim + 2 ] += fz / 3.0;

        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void FPP_4 ( double *F, int **f_pre_ind, double *f_pre_val, int FPP,
             double **coord, int dim )
{
    /** Descripción **/
    // Impone las fuerzas de presión al vector de fuerzas. Si es positiva
    // empuja la cara y si es negativa la jala. Esto para elementos tipo 4
    // (hexaedro) en 3D.

    /** Variables de entrada **/
    // {F}          Vector de fuerzas
    // FP           Número de fzas. de presión.
    // {f_pre_ind}  Caras en las que está aplicada la fuerzas.
    // {f_punt_val} Valores de las fuerzas de presión.
    // [coord]      Coordenadas de todos los nodos.
    // dim          Dimensión.

    /** Variables locales **/
    int i, a, b, c, j, d;
    double ax, ay, az, bx, by, bz, cx, cy, cz, a_b, b_c, Ma, Mb, Mc;
    double Pab, Pcb, h1, h2, area, nx, ny, nz, Mn, f, fx, fy, fz;

    for ( i = 0; i < FPP; i++ ){

        a = f_pre_ind [ i ][ 0 ];
        b = f_pre_ind [ i ][ 1 ];
        c = f_pre_ind [ i ][ 2 ];
        d = f_pre_ind [ i ][ 3 ];

        ax = coord[ b ][ 0 ] - coord[ a ][ 0 ];
        ay = coord[ b ][ 1 ] - coord[ a ][ 1 ];
        az = coord[ b ][ 2 ] - coord[ a ][ 2 ];

        bx = coord[ c ][ 0 ] - coord[ a ][ 0 ];
        by = coord[ c ][ 1 ] - coord[ a ][ 1 ];
        bz = coord[ c ][ 2 ] - coord[ a ][ 2 ];

        cx = coord[ d ][ 0 ] - coord[ a ][ 0 ];
        cy = coord[ d ][ 1 ] - coord[ a ][ 1 ];
        cz = coord[ d ][ 2 ] - coord[ a ][ 2 ];

        //Producto punto a.b y b.c
        a_b = ax * bx + ay * by + az * bz;
        b_c = bx * cx + by * cy + bz * cz;

        // Magnitud de a, b y c
        Ma = sqrt ( ax * ax + ay * ay + az * az );
        Mb = sqrt ( bx * bx + by * by + bz * bz );
        Mc = sqrt ( cx * cx + cy * cy + cz * cz );

        // Proyección de a sobre b
        Pab = a_b / Mb;

        // Proyección de a sobre c sobre b
        Pcb = b_c / Mb;

        // Alturas
        h1 = sqrt ( Ma * Ma - Pab * Pab );
        h2 = sqrt ( Mc * Mc - Pcb * Pcb );

        // Area del triángulo
        area = 0.5 * Mb * ( h1 + h2 );

        // Fuerza total
        f = f_pre_val[ i ] * area;

        // Producto cruz a x b
        nx = ay * bz - az * by;
        ny = az * bx - ax * bz;
        nz = ax * by - ay * bx;

        // Magnitud del vector normal
        Mn = sqrt ( nx * nx + ny * ny + nz * nz );

        // Vectores unitarios en x, y, z
        nx /= Mn;
        ny /= Mn;
        nz /= Mn;

        // Componentes del vector f
        fx = f * nx;
        fy = f * ny;
        fz = f * nz;

        // Suma de componente al vector F
        for ( j = 0; j < 4; j++ ){

            a = f_pre_ind[ i ][ j ];

            F[ a * dim + 0 ] += fx / 4.0;
            F[ a * dim + 1 ] += fy / 4.0;
            F[ a * dim + 2 ] += fz / 4.0;

        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void neu_FS ( double *F, int **f_sup_ind, double **f_sup_val, int FS,
              double **coord, int tipo_elem, int nodos_elem, int dim )
{
    /** Descripción **/
    // Impone las fuerzas de superficie al vector de fuerzas.

    /** Variables de entrada **/
    // {F}          Vector de fuerzas
    // FS           Número de fuerzas de superficie.
    // {f_sup_ind}  Caras en las que están plicadas las fuerzas.
    // {f_sup_val}  Valores de las fuerzas
    // tipo_elem    Tipo de elemento
    // [coord]      Coordenadas de todos los nodos
    // dim          Dimensión
    // tipo_elem    Tipo de elemento de la malla


    /** Variables locales **/

    if ( tipo_elem == 3 ) FS_3 ( F, f_sup_ind, f_sup_val, FS, coord,
                                 dim );

    if ( tipo_elem == 4 ) FS_4 ( F, f_sup_ind, f_sup_val, FS, coord,
                                 dim );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void FS_3 ( double *F, int **f_sup_ind, double **f_sup_val, int FS,
            double **coord, int dim )
{
    /** Descripción **/
    // Impone las fuerzas de superficie al vector de fuerzas. Esto para
    // elementos tipo 3 (tetraedros) en 3D.

    /** Variables de entrada **/
    // {F}          Vector de fuerzas
    // FS           Número de fuerzas de superficie.
    // {f_sup_ind}  Caras en las que están plicadas las fuerzas.
    // {f_sup_val}  Valores de las fuerzas
    // [coord]      Coordenadas de todos los nodos
    // dim          Dimensión

    /** Variables locales **/
    int i, a, b, c, j;
    double ax, ay, az, bx, by, bz, a_b, Ma, Mb, Pab, h, area;
    double fx, fy, fz;

    for ( i = 0; i < FS; i++ ){

        a = f_sup_ind [ i ][ 0 ];
        b = f_sup_ind [ i ][ 1 ];
        c = f_sup_ind [ i ][ 2 ];

        ax = coord[ b ][ 0 ] - coord[ a ][ 0 ];
        ay = coord[ b ][ 1 ] - coord[ a ][ 1 ];
        az = coord[ b ][ 2 ] - coord[ a ][ 2 ];

        bx = coord[ c ][ 0 ] - coord[ a ][ 0 ];
        by = coord[ c ][ 1 ] - coord[ a ][ 1 ];
        bz = coord[ c ][ 2 ] - coord[ a ][ 2 ];

        //Producto punto a.b
        a_b = ax * bx + ay * by + az * bz;

        // Magnitud de a y b
        Ma = sqrt ( ax * ax + ay * ay + az * az );
        Mb = sqrt ( bx * bx + by * by + bz * bz );

        // Proyección de a sobre b
        Pab = a_b / Mb;

        // Altura
        h = sqrt ( Ma * Ma - Pab * Pab );

        // Area del triángulo
        area = 0.5 * Mb * h;

        // Componentes de la fuerza
        fx = f_sup_val[ i ][ 0 ] * area;
        fy = f_sup_val[ i ][ 1 ] * area;
        fz = f_sup_val[ i ][ 2 ] * area;

        // Suma de componente al vector F
        for ( j = 0; j < 3; j++ ){

            a = f_sup_ind[ i ][ j ];

            F[ a * dim + 0 ] += fx / 3.0;
            F[ a * dim + 1 ] += fy / 3.0;
            F[ a * dim + 2 ] += fz / 3.0;

        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void FS_4 ( double *F, int **f_sup_ind, double **f_sup_val, int FS,
            double **coord, int dim )
{
    /** Descripción **/
    // Impone las fuerzas de superficie al vector de fuerzas. Esto para
    // elementos tipo 4 (hexaedro) en 3D.

    /** Variables de entrada **/
    // {F}          Vector de fuerzas
    // FS           Número de fuerzas de superficie.
    // {f_sup_ind}  Caras en las que están plicadas las fuerzas.
    // {f_sup_val}  Valores de las fuerzas
    // [coord]      Coordenadas de todos los nodos
    // dim          Dimensión

    /** Variables locales **/
    int i, a, b, c, j, d;
    double ax, ay, az, bx, by, bz, cx, cy, cz, a_b, b_c, Ma, Mb, Mc;
    double Pab, Pcb, h1, h2, area, fx, fy, fz;

    for ( i = 0; i < FS; i++ ){

        a = f_sup_ind [ i ][ 0 ];
        b = f_sup_ind [ i ][ 1 ];
        c = f_sup_ind [ i ][ 2 ];
        d = f_sup_ind [ i ][ 3 ];

        ax = coord[ b ][ 0 ] - coord[ a ][ 0 ];
        ay = coord[ b ][ 1 ] - coord[ a ][ 1 ];
        az = coord[ b ][ 2 ] - coord[ a ][ 2 ];

        bx = coord[ c ][ 0 ] - coord[ a ][ 0 ];
        by = coord[ c ][ 1 ] - coord[ a ][ 1 ];
        bz = coord[ c ][ 2 ] - coord[ a ][ 2 ];

        cx = coord[ d ][ 0 ] - coord[ a ][ 0 ];
        cy = coord[ d ][ 1 ] - coord[ a ][ 1 ];
        cz = coord[ d ][ 2 ] - coord[ a ][ 2 ];

        //Producto punto a.b y b.c
        a_b = ax * bx + ay * by + az * bz;
        b_c = bx * cx + by * cy + bz * cz;

        // Magnitud de a, b y c
        Ma = sqrt ( ax * ax + ay * ay + az * az );
        Mb = sqrt ( bx * bx + by * by + bz * bz );
        Mc = sqrt ( cx * cx + cy * cy + cz * cz );

        // Proyección de a sobre b
        Pab = a_b / Mb;

        // Proyección de a sobre c sobre b
        Pcb = b_c / Mb;

        // Alturas
        h1 = sqrt ( Ma * Ma - Pab * Pab );
        h2 = sqrt ( Mc * Mc - Pcb * Pcb );

        // Area del triángulo
        area = 0.5 * Mb * ( h1 + h2 );

        // Componentes de la fuerza
        fx = f_sup_val[ i ][ 0 ] * area;
        fy = f_sup_val[ i ][ 1 ] * area;
        fz = f_sup_val[ i ][ 2 ] * area;

        // Suma de componente al vector F
        for ( j = 0; j < 4; j++ ){

            a = f_sup_ind[ i ][ j ];

            F[ a * dim + 0 ] += fx / 4.0;
            F[ a * dim + 1 ] += fy / 4.0;
            F[ a * dim + 2 ] += fz / 4.0;

        }

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
