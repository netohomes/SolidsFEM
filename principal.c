//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

        /** MECANICA DE SOLIDOS RESUELTO POR ELEMENTOS FINITOS **/
                             // Para 3D

    // Elementos 3D: tetraedros(3) y hexaedros(4).
    // Dim. No.      Tipo de elemento    P.G.    C.C v�lidas
    //                                           Fuerzas     Despl.
    // 3    4        Tetraedro           1       P,L,S,V     P,L,S
    // 3    5        Hexaedro            8       P,L,S,V     P,L,S
    // Nota: P -> Sobre puntos, L -> Sobre lineas, S -> Sobre superficies
    //       V -> Sobre volumenes
    // El numero de nodos no debe llegar a 100 millones

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "principal.h"
#include "lectura.h"
#include "memoria.h"
#include "matriz_global.h"
#include "cond_cont.h"
#include "grad_conj.h"
#include "result.h"
#include "esfuerzos.h"

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

int main ( int argc, const char *argv[] )
{
    //====================================================================
     /** Descripcion **/
        // Funcion principal del programa
    //====================================================================
    /** Variables de entrada **/
        // argc     Numero de argumentos para el programa
        // argv[]   Direccion al archivo con el que se esta trabajando
    //====================================================================
        lectura ( (char *)argv[ 1 ] );
    //====================================================================
        ady ( elemento, n_elem, nodos, nodos_elem );
    //====================================================================
        rigid ( tabla_ady, elemento, f_cpo_ind, f_cpo_val, tipo_elem,
                nodos_elem, pg_elem, dim );
    //====================================================================
        cond_cont ( K, F, tabla_ady, coord, nodos, tipo_elem, nodos_elem,
                    FP, f_punt_ind, f_punt_val, FPP, f_pre_ind, f_pre_val,
                    FS, f_sup_ind, f_sup_val, DL, d_lin_ind, d_lin_val,
                    DS, d_sup_ind, d_sup_val, dim );
    //====================================================================
        grad_conj ( K, F, tabla_ady, U, nodos * dim, itera, error, dim );
    //====================================================================
        esfuerzos ( U, coord, elemento, dim, nodos, n_elem, tipo_elem,
                    nodos_elem );
    //====================================================================
        result ((char *)argv[ 1 ], n_elem, pg_elem, tipo_elem, nodos, U,
                ESF, dim);
    //====================================================================

    return 0;

}

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

