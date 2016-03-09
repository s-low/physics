/*
 ============================================================================
 Name        : Diffusion.c
 Author      : Sam Low
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

int main(void) {

	double lambda, alpha, delta_t, h, L, t, t_max;
	int i, j, n, s;
	FILE *datafile;

	printf("How many seconds to iterate:");
	scanf("%lf", &t_max);

	printf("\nPoker is 50cm long\n");
	printf("How many nodes along the poker:");
	scanf("%d", &n);

	L = 0.5; //length in metres
	delta_t = 1; //time increment = 1 second
	alpha = 1.66 * pow(10, -5); //diffusion coefficient

	h = L / n; //distance between nodes

	gsl_matrix * A = gsl_matrix_alloc(n, n);
	gsl_vector * phi_0 = gsl_vector_alloc(n);
	gsl_vector * phi_t = gsl_vector_alloc(n);
	gsl_permutation * p = gsl_permutation_alloc(n); //permutation matrix. required for LU decomp.
	// A(mat) * phi_t(vec) = phi_0(vec)

	lambda = (alpha * delta_t) / pow(h, 2); //constant lambda. see report.

	//setting up matrix A
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {

			if (i != 0 && i != n - 1) {
				gsl_matrix_set(A, i, i, 1 + (2 * lambda));
				gsl_matrix_set(A, i, i - 1, -lambda);
				gsl_matrix_set(A, i, i + 1, -lambda);
			}

			else if (i == 0) {
				gsl_matrix_set(A, i, i, 1 + lambda);
				gsl_matrix_set(A, i, i + 1, -lambda);
			}

			else if (i == n - 1) {
				gsl_matrix_set(A, i, i, 1 + lambda);
				gsl_matrix_set(A, i, i - 1, -lambda);
			}

			else {
				gsl_matrix_set(A, i, j, 0); //all other elements zero
			}
		}

	//poker initially at 20C
	for (i = 0; i < n; i++) {
		gsl_vector_set(phi_0, i, 20);
		gsl_vector_set(phi_t, i, 0);
	}

	gsl_vector_set(phi_0, n - 1, 1000); //poker tip in furnace

	gsl_linalg_LU_decomp(A, p, &s);

	datafile = fopen("Temp.txt", "w");
	if (datafile != NULL) //check to see if the file was actually opened
	{
		for (t = 0; t <= t_max; t += delta_t) { //looping forwards in time

			gsl_linalg_LU_solve(A, p, phi_0, phi_t); //solve for new phi's. now in phi_t.

			gsl_vector_set(phi_0, n - 1, 1000);
			gsl_vector_set(phi_t, n - 1, 1000); //reset poker tip

			for (i = 0; i < n; i++) {
				gsl_vector_set(phi_0, i, gsl_vector_get(phi_t, i));
			} //copy phi_t back to phi_0. redo the process.
		}

		for (i = 0; i < n; i++) {
			fprintf(datafile, "%f\t%lf\n", i * h, gsl_vector_get(phi_t, i));
		}
		fprintf(datafile, "\n");

	}
}
