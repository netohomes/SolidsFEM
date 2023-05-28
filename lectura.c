//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                      /** LECTURA DE DATOS **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "principal.h"
#include "lectura.h"
#include "memoria.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

void lectura ( char *r )
{

    //====================================================================
    /** Descripcion **/
    // Lee los datos proporcionados por GID. Estos son:
    //  - Coordenadas de los nodos
    //  - Conectividades de los elementos
    //  - Propiedades de los materiales
    //  - Condiciones de contorno
    //====================================================================
    /** Variables de entrada **/
    // r            Contiene la direccion del archivo de datos
    //              (sin extension)
    // {elemento}   Vector que guarda las propiedades de los elementos y
    //              sus conectividades
    // [coord]      Matriz de coordenadas
    //              Variables auxuliares para guardar:
    // dima_a       La dimension
    // tipo_elem_a  El tipo de elemento empleado para mallar
    // n_elem_a     N�mero de elementos
    // nodos_a      Nodos del elemento
    //====================================================================
    /** Variables locales **/
    // ruta Contendra la ruta completa del archivo de datos
    // aux  Empleada para guardar un valor que no es empleado a lo largo
    //      del programa

    char *ruta;
    int i, j, aux;
    FILE *fp;

    //====================================================================

    ruta = ( char *) calloc ( 256, sizeof ( char ) );
    strcpy ( ruta, r );
    strcat ( ruta, ".dat" );

    fp = fopen ( ruta , "r" );
	if (fp == NULL){
		printf ( "Archivo de datos no encontrado\n" );
		exit ( 1 );
	}

	//====================================================================

    fscanf ( fp, "===================================================" );
    fscanf ( fp, "================\n" );
    fscanf ( fp, "                         ARCHIVO DE DATOS\n" );
    fscanf ( fp, "===================================================" );
    fscanf ( fp, "================\n" );

    //====================================================================

    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );

    fscanf ( fp, "DIMENSION   %d\n", &dim       );
    fscanf ( fp, "TIPO_ELEM   %d\n", &tipo_elem );
    fscanf ( fp, "ELEMENTOS   %d\n", &n_elem    );
    fscanf ( fp, "NODOS       %d\n", &nodos     );

    //====================================================================

    tipo_elem -= 1;
    nodos_elem = prop_elemento ( tipo_elem, dim, 1 );
    pg_elem    = prop_elemento ( tipo_elem, dim, 2 );

    //====================================================================

    elemento = mem_vector_elem ( n_elem, "elementos." );
    coord    = mem_matriz ( nodos, dim, "coordendas." );

    for ( i = 0; i < n_elem; i++ ){
        elemento[i].conec = mem_vector_int (nodos_elem,"de conectividad.");
	}

    //====================================================================

    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "CONECTIVIDAD\n" );
    if ( tipo_elem == 4 ){
        fscanf ( fp, "    ELEM.         NODOS                      " );
        fscanf ( fp, "                                             " );
        fscanf ( fp, "         E(kg/cm2)  v\n" );
    }
    if ( tipo_elem == 3 ){
        fscanf ( fp, "    ELEM.         NODOS                      " );
        fscanf ( fp, "   E(kg/cm2)  v\n" );
    }

    for ( i = 0; i < n_elem; i++ ){

        fscanf ( fp, "%d", &aux );

        for ( j = 0; j < nodos_elem; j++ ){
            fscanf ( fp, "%d", &elemento[ i ].conec[ j ] );
            elemento[ i ].conec[ j ] -= 1;
        }

        fscanf ( fp, "%lf", &elemento[ i ].E );
        fscanf ( fp, "%lf\n", &elemento[ i ].v );

    }

    //====================================================================

    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "COORDENADAS\n" );
    fscanf ( fp, "NODO                    X              Y          " );
    fscanf ( fp, "    Z\n" );

    for ( i = 0; i < nodos; i++ ){

        fscanf ( fp, "%d", &aux );

        if ( dim == 3 ) aux = 3;
        for ( j = 0; j < aux; j++ ){
            fscanf ( fp, "%lf", &coord[ i ][ j ] );
            if ( j == aux - 1 ) fscanf ( fp, "\n" );
        }

    }

    //====================================================================

    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "FUERZAS PUNTUALES\n" );
    fscanf ( fp, "CANTIDAD   %d\n", &FP );
    fscanf ( fp, "NODO                   FX              FY          " );
    fscanf ( fp, "    FZ\n");

    f_punt_ind = mem_vector_int ( FP, "indices de fzas. puntuales." );
    f_punt_val = mem_matriz     ( FP, dim, "valores de fzas puntuales." );

    for ( i = 0; i < FP; i++ ){

        fscanf ( fp, "%d", &f_punt_ind[ i ] );
        f_punt_ind[ i ] -= 1;

        for ( j = 0; j < dim; j++ ){

            fscanf ( fp, "%lf", &f_punt_val[ i ][ j ] );

        }

        fscanf ( fp, "\n" );

    }

    //====================================================================

    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "FUERZAS DE PRESION\n" );
    fscanf ( fp, "CANTIDAD   %d\n", &FPP );

    if ( tipo_elem == 3 ){
        fscanf ( fp, "NODOS                                     F\n"     );
        aux = 3;
    }

    if ( tipo_elem == 4 ){
        fscanf ( fp, "NODOS                                          " );
        fscanf ( fp, "    F\n" );
        aux = 4;
    }

    f_pre_ind = mem_matriz_int ( FPP, aux, "indices de f_pre." );
    f_pre_val = mem_vector     ( FPP, "valores de f_pre." );

    for ( i = 0; i < FPP; i++ ){

        for ( j = 0; j < aux; j++ ){

            fscanf ( fp, "%d", &f_pre_ind[ i ][ j ] );
            f_pre_ind[ i ][ j ] -= 1;

        }

        fscanf ( fp, "%lf\n", &f_pre_val[ i ] );

    }

    //====================================================================

    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "FUERZAS DE SUPERFICIE\n" );
    fscanf ( fp, "CANTIDAD   %d\n", &FS );

    if ( tipo_elem == 3 ){
        fscanf ( fp, "NODOS                                    FX    " );
        fscanf ( fp, "          FY              FZ\n" );
        aux = 3;
    }

    if ( tipo_elem == 4 ){
        fscanf ( fp, "NODOS                                          " );
        fscanf ( fp, "  FX               FY              FZ\n" );
        aux = 4;
    }

    f_sup_ind = mem_matriz_int ( FS, aux, "indices de f_sup." );
    f_sup_val = mem_matriz     ( FS, dim, "valores de f_sup." );

    for ( i = 0; i < FS; i++ ){

        for ( j = 0; j < aux; j++ ){

            fscanf ( fp, "%d", &f_sup_ind[ i ][ j ] );
            f_sup_ind[ i ][ j ] -= 1;

        }

        for ( j = 0; j < dim; j++ ){

            fscanf ( fp, "%lf", &f_sup_val[ i ][ j ] );

        }

        fscanf ( fp, "\n" );

    }

    //====================================================================

    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "FUERZAS DE CUERPO\n" );
    fscanf ( fp, "CANTIDAD   %d\n", &FC );

    f_cpo_ind = mem_vector_int ( FC, "indices de fzas. de cuerpo." );
    f_cpo_val = mem_matriz     ( n_elem, dim, "valores de fzas. de cuerpo." );
    inicializa_matriz ( f_cpo_val, n_elem, dim );

    fscanf ( fp, "ELEM.                  FX              FY           " );
    fscanf ( fp, "   FZ\n" );

    for ( i = 0; i < FC; i++ ){

        fscanf ( fp, "%d", &f_cpo_ind[ i ] );
        f_cpo_ind[ i ] -= 1;

        for ( j = 0; j < dim; j++ ){

            fscanf ( fp, "%lf", &f_cpo_val[ i ][ j ] );

        }

        fscanf( fp, "\n" );

    }

    //====================================================================
    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "DESPLAZAMIENTOS EN LINEA\n" );
    fscanf ( fp, "CANTIDAD   %d\n", &DL );

    fscanf ( fp, "NODOS                  DX              DY           " );
    fscanf ( fp, "   DZ\n" );

    d_lin_ind = mem_vector_int ( DL, "indices de despl. en lineas." );
    d_lin_val = mem_matriz     ( DL, dim, "valores de despl. en lineas" );

    for ( i = 0; i < DL; i++ ){

        fscanf ( fp, "%d", &d_lin_ind[ i ] );
        d_lin_ind[ i ] -= 1;

        for ( j = 0; j < dim; j++ ){

            fscanf ( fp, "%lf", &d_lin_val[ i ][ j ] );

        }

        fscanf( fp, "\n" );

    }

    //====================================================================
    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );
    fscanf ( fp, "DESPLAZAMIENTOS EN SUP.\n" );
    fscanf ( fp, "CANTIDAD   %d\n", &DS );

    if ( tipo_elem == 3 ){
        fscanf ( fp, "NODOS                                    DX    " );
        fscanf ( fp, "          DY              DZ\n" );
        aux = 3;
    }

    if ( tipo_elem == 4 ){
        fscanf ( fp, "NODOS                                          " );
        fscanf ( fp, "  DX               DY              DZ\n" );
        aux = 4;
    }

    d_sup_ind = mem_matriz_int ( DS, aux, "indices de despl. de sup." );
    d_sup_val = mem_matriz     ( DS, dim, "valores de despl. de sup." );

    for ( i = 0; i < DS; i++ ){

        for ( j = 0; j < aux; j++ ){

            fscanf ( fp, "%d", &d_sup_ind[ i ][ j ] );
            d_sup_ind[ i ][ j ] -= 1;

        }

        for ( j = 0; j < dim; j++ ){

            fscanf ( fp, "%lf", &d_sup_val[ i ][ j ] );

        }

        fscanf ( fp, "\n" );

    }
    //====================================================================
    fscanf ( fp, "..................................................." );
    fscanf ( fp, "................\n" );

    fscanf ( fp, "ITERACIONES %d",  &itera );
    fscanf ( fp, "ERROR       %lf", &error );

    //====================================================================

    fclose ( fp );

}

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

int prop_elemento ( int tipo_elem, int dim, int prop )
{
    /** Descripci�n **/
    // Determina las caracter�sticas del elemento que se us� para mallar

    int nodos_elem, pg_elem;

    if ( dim == 3 ){

        switch ( tipo_elem ){

            case 3:
            {
                nodos_elem = 4;
                pg_elem    = 1;
                break;
            }

            case 4:
            {
                nodos_elem = 8;
                pg_elem    = 8;
                break;
            }

            default:
            {
                printf ( "\nError: no se pudo determinar nodos_elem y pg_elem. " );
                printf ( "\nCaso no incluido para 3D." );
                exit ( 1 );
            }

        }

    }

    if ( prop == 1 ){
        return nodos_elem;
    }
    else{
        return pg_elem;
    }

}

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}


