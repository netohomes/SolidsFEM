//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

                     /**  OPERACIONES MATRICIALES Y VECTORIALES  **/

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

#include "oper_matrices.h"
#include <stdio.h>
#include <stdlib.h>

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void multiplica ( double **C, double **A, double **B, int rengA, int col, int colB )
{
	/** Descripci�n **/
	// Multiplica [A] x [B] y el resultado lo guarda en [C]

	int i, j, k;
	for ( i = 0; i < rengA; i++ ){
		for ( j = 0; j < colB; j++ ){
			for ( k = 0; k < col; k++ ){
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void mult_mat_vector ( double *C, double **A, double *B, int rengA, int colA, int rengB )
{
	/** Descripci�n **/
	// Realiza la multiplicaci�n matriz - vector. Multiplica [A] x {B} y el resultado lo guarda
	// en {C}

	int i, j;

	if ( colA != rengB ) puts ( "Error: no se puede realizar multiplicacion mat X vector." );

	for ( i = 0; i < rengA; i++ ){
		for ( j = 0; j < rengB; j++ ){
			C[i] += A[i][j]*B[j];
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void transponer_mat ( double **A, double **B, int reng, int col )
{
    /** Descripci�n **/
	// Transpone una matriz [A] y el resuldato lo guarda en otra matriz [B]
	// reng y col son los renglones y las columnas de A.

	int i, j;
	for ( i = 0; i < reng; i++ ){
		for ( j = 0; j < col; j++ ){
			B[j][i] = A[i][j];
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void mult_mat_escalar ( double **A, double b, int reng, int col)
{
	/** Descripci�n **/
	// Multiplica una matriz [A] por un escalar b y el resultado es guardado en la misma matriz [A]

	int i, j;
	for ( i = 0; i < reng; i++ ){
		for ( j = 0; j < col; j++){
			A[i][j] *= b;
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void mult_vec_escalar ( double *v, double b, int reng )
{
    /** Descripci�n **/
	// Multiplica un vector {v} por un escalar b y el resultado es guardado en el mismo vector {v}

	int i;
	for ( i = 0; i < reng; i++ ){
        v[i] *= b;
	}

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void suma_mat ( double **A, double **B, int reng, int col)
{
	/** Descripci�n **/
	// Suma dos matrices [A] + [B] y el resultado los guarda en [A]

	int i, j;
	for ( i = 0; i < reng; i++ ){
		for ( j = 0; j < col; j++){
			A[i][j] += B[i][j];
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void suma_vec ( double *a, double *b, int reng )
{
	/** Descripci�n **/
	// Suma dos vectores {a} + {b} y el resultado los guarda en {a}

	int i;
	for ( i = 0; i < reng; i++ ){
		a[ i ] += b[ i ];
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
