//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

        /** SOLVER: GRADIENTE CONJUGADO PARA MATRICES RALAS **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "principal.h"
#include "grad_conj.h"
#include "memoria.h"
#include <stdio.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void grad_conj ( double **K, double *F, struct adyacencias *A, double *x,
                 int n, int n_iter, double tol, int dim )
{

    /** Descripción **/
    // Resuelve un sistema de ecuaciones con la matriz K guardada en formato
    // de renglón comprimido

    /** Variables de entrada **/
    // [K]      Matriz en formato de renglón comprimido
    //          Su numeración tiene que empezar en 0.
    // {F}      Vector de fuerzas.
    // n        Tamaño original de la matriz ( 3 x nodos )
    // tol      Tolerancia
    // n_iter   Cantidad máxima de iteraciones a realizar

    /** Variables de salida **/
    // x        Vector solución

    /** Variables locales **/
    int i, j, k, a, b;
    double *g, *p, *w;
    double epsilon, alfa, beta;
    double suma, suma_g1, suma_g2;
    suma_g1 = 0.0;
    suma_g2 = 0.0;

    g = mem_vector ( n, "g");
    p = mem_vector ( n, "p");
    w = mem_vector ( n, "w");

    for ( i = 0; i < n; i++ ){
        g[ i ] = 0.0;
        p[ i ] = 0.0;
        w[ i ] = 0.0;
    }

    epsilon = tol * tol;

    // Vector inicial: {x0}
    for ( i = 0; i < n; i++) x[ i ] = 0.0;

    // Gradiente inicial: {p0} = {-g0} = [A]{x0} - b
    for ( i = 0; i < n; i++ ){
        suma = 0.0;
        a = ( i - i % dim ) / dim;
        for ( j = 0; j < A[ a ].cols * dim; j++ ){
            b = ( j - j % dim ) / dim;
            b = A[ a ].ind[ b ];
            b = b * dim + j % dim;
            suma +=  K[ i ][ j ] * x[ b ];
        }
        g[ i ] = suma - F[ i ];//Posible error
        p[ i ] = - g[ i ];
        suma_g1 += g[ i ] * g[ i ];
    }

    k = 0;
    while ( k < n_iter )
    {

        if ( suma_g1 <= epsilon) break;

        // {w} = [A]{pk}
        for ( i = 0; i < n; i++ ){
            suma = 0;
             a = ( i - i % dim ) / dim;
            for ( j = 0; j < A[ a ].cols * dim; j++ ){
                b = ( j - j % dim ) / dim;
                b = A[ a ].ind[ b ];
                b = b * dim + j % dim;
                suma += K[ i ][ j ] * p[ b ];
            }
            w[ i ] = suma;
        }

        // alfa_k = {gk_T}{gk}/({pk_T}{w})
        suma = 0.0;
        for ( i = 0; i < n; i++ ){
            suma += p[ i ] * w[ i ];
        }
        alfa = suma_g1 / suma;

        // {xk+1} = {xk} + alfa_k{pk}
        for ( i = 0; i < n; i++){
            x[ i ] = x[ i ] + alfa * p[ i ];// x redefinido
            g[ i ] = g[ i ] + alfa * w[ i ];// g redefinido
        }

        // beta_k = {gk+1_T}{gk}/({gk_T}{gk})
        suma_g2 = 0.0;
        for ( i = 0; i < n; i++ ){
            suma_g2 += g[ i ] * g[ i ];
        }
        beta = suma_g2 / suma_g1;

        // {pk+1} = -{gk+1}+beta_k{pk}
        for ( i = 0; i < n; i++ ){
            p[ i ] = -g[ i ] + beta * p[ i ]; // p redefinido
        }

        suma_g1 = suma_g2;

        k++;
    }

    libera_vector ( g );
    libera_vector ( p );
    libera_vector ( w );

    /*printf ( "\n\nVector solucion" );
    for ( i = 0; i < n; i++ ){
        printf ( "\n%.6lf", x[ i ] );
    }*/

}

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
