/*
 ============================================================================
 Name        : SVD.c
 Author      : Sam Low

 ============================================================================
 */

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>

/*Program will request dimension n, then print an nxn matrix of
 random numbers to screen, followed by its inverse */

int main(void) {
	int n; //dimension of matrix
	int i, j;

	printf("Enter the number of rows and columns n:");
	scanf("%d", &n);

	gsl_matrix * A = gsl_matrix_alloc(n, n);
	gsl_matrix * V = gsl_matrix_alloc(n, n);
	gsl_matrix * Sm = gsl_matrix_alloc(n, n);
	gsl_matrix * X = gsl_matrix_alloc(n, n);
	gsl_matrix * Y = gsl_matrix_alloc(n, n);

	gsl_vector * S = gsl_vector_alloc(n); //only the diagonal elements of S
	gsl_vector * work = gsl_vector_alloc(n); //workspace for GSL decomp

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			gsl_matrix_set(A, i, j, (20 * rand()) / (double) RAND_MAX);
		}

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			printf("A(%d,%d) = %g\n", i, j, gsl_matrix_get(A, i, j));

	printf("\n\n");

	gsl_linalg_SV_decomp(A, V, S, work);

	// A decomposed into: A = U S V
	// A now replaced by U
	// V NOT TRANSPOSED by GSL yet
	// find Pseudo inverse of S, then calculate: A(PI) = V S(PI) U^T

	for (i = 0; i < n; i++) {
		if (gsl_vector_get(S, i) != 0) {
			gsl_matrix_set(Sm, i, i, (1 / gsl_vector_get(S, i)));
		} else {
			gsl_matrix_set(Sm, i, i, (gsl_vector_get(S, i)));
		}
	}

//now Sm is matrix with 1/S elements on diagonal
//V * S * U^T is now the aim

//FIND: A(Ps.Inv) = V Sm(Ps.Inv) U^T
// X = Sm*U
// Y = V*X = V*Sm*U

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Sm, A, 0.0, X); //multiply Sm.A, put in X
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, X, 0.0, Y); //multiply V.X, put in Y


	//Y now stores the inverse
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("Inv(%d,%d) = %g\n", i, j, gsl_matrix_get(Y, i, j));
		}
	}

}
