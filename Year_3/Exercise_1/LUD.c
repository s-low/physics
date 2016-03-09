/*
 ============================================================================
 Name        : LUD.c
 Author      : Sam Low

 ============================================================================
 */

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

/*Program will request dimension n, then print an nxn matrix of
 random numbers to screen, followed by its inverse */

int main(void) {
	int n;
	int s; //signum s for LU decomp
	int i, j;

	printf("Enter the number of rows and columns n:");
	scanf("%d", &n);

	gsl_matrix * A = gsl_matrix_alloc(n, n);
	gsl_matrix * inverse = gsl_matrix_alloc(n, n);
	gsl_permutation * p = gsl_permutation_alloc(n); //perm.matrix needed to reorder

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			gsl_matrix_set(A, i, j, (20 * rand()) / (double) RAND_MAX);

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			printf("A(%d,%d) = %g\n", i, j, gsl_matrix_get(A, i, j));

	gsl_linalg_LU_decomp(A, p, &s);
	gsl_linalg_LU_invert(A, p, inverse);

	printf("\n\n");

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			printf("Inv(%d,%d) = %g\n", i, j, gsl_matrix_get(inverse, i, j));

}
