/*
 ============================================================================
 Name        : Jacobi.c
 Author      : Sam Low
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*edges is used to set all edges of matrix to same value
 * allocate is used to allocate a 2d array*/

void edges(int dim1, int dim2, double ***m, double value);
void allocate(int dim_x, int dim_y, double ***m);

int main() {

	int i, i1 = 0, j, loop;
	double **V, **Vnew; //jacobi method - two complete copies of grid
	int dimx, dimy; //grid dimensions
	double diff, perc, X;
	FILE *datafile;

	printf("Program uses plate capacitor with extremal corner of y=60cm, x=60cm\n");
	printf("Be sure to exceed these dimensions\n");

	printf("\nx dimension in whole centimetres:");
	scanf("%d", &dimx);

	printf("\ny dimension in whole centimetres:");
	scanf("%d", &dimy);

	printf("\nEnter the the convergence condition in percent:");
	scanf("%lf", &X);
	//convergence criterion. Stop iterating if no node change by more than X%

	allocate(dimx, dimy, &V);
	allocate(dimx, dimy, &Vnew);

	for (i = 0; i < dimy; i++) { //fill V with an initial guess
		for (j = 0; j < dimx; j++) {
			V[i][j] = 0;
		}
	}

	/*BOUNDARY CONDITIONS: In this case, a plate capacitor at y=40,60 x=40-60*/

	edges(dimx, dimy, &V, 0); //set all edges to zero

	for (i = 40; i < 60; i++) { //+ve plate
		V[i][60] = 10;
	}

	for (i = 40; i < 60; i++) { //-ve plate
		V[i][40] = -10;
	}

	//now iterate finite difference method

	datafile = fopen("Jacobi.txt", "w");
	if (datafile != NULL) //check to see if the file was actually opened
	{

		loop = 1;
		while (loop == 1) { //implementing the convergence criterion as a while loop

			/*everything here down iterates over the whole array (once)
			 * Is enclosed in while loop that keeps iterating if a single
			 * element changes by more than X%*/

			loop = 0;
			for (i = 1; i < dimy - 1; i++)
				for (j = 1; j < dimx - 1; j++) {
					Vnew[i][j] = 0.25
							* ((V[i + 1][j]) + (V[i - 1][j]) + (V[i][j + 1])
									+ (V[i][j - 1]));
				}

			for (i = 40; i < 60; i++) { //+ve plate
				Vnew[i][60] = 10;
				V[i][60] = 10;
			}

			for (i = 40; i < 60; i++) { //-ve plate
				Vnew[i][40] = -10;
				V[i][40] = -10;
			}

			edges(dimx, dimy, &V, 0);
			edges(dimx, dimy, &Vnew, 0);

			for (i = 1; i < dimy - 1; i++)
				for (j = 1; j < dimx - 1; j++) {
					diff = pow(V[i][j] - Vnew[i][j], 2);
					perc = 100 * (sqrt(diff) / ((V[i][j] + Vnew[i][j]) / 2));

					if (perc > X) { //convergence criterion: no node changes > 0.1%
						loop = 1;
					}
				}

			for (i = 0; i < dimy; i++) { //overwrite V with new values
				for (j = 0; j < dimx; j++) {
					V[i][j] = Vnew[i][j];
				}
			}

		}

		for (i = 0; i < dimy; i++)  //overwrite V with new values
			for (j = 0; j < dimx; j++) {

				if (i != i1) {
					fprintf(datafile, "\n"); //quirk of using gnuplot with pm3d
				}
				i1 = i;

				fprintf(datafile, "%d\t%d\t%lf\n", i, j, V[i][j]);

			}

		fclose(datafile);
	}
}

void edges(int dim_x, int dim_y, double ***m, double value) {

	int i, j;

	for (i = 0; i < dim_y; i++) { //set right edge to argument value
		(*m)[i][dim_x - 1] = value;
	}

	for (j = 0; j < dim_x; j++) { //set top edge
		(*m)[0][j] = value;
	}

	for (j = 0; j < dim_x; j++) { //set bottom edge
		(*m)[dim_y - 1][j] = value;
	}

	for (i = 0; i < dim_y; i++) { //set left edge
		(*m)[i][0] = value;
	}

}

void allocate(int dim_x, int dim_y, double ***m) {
	int i;

	(*m) = malloc(dim_y * sizeof(double*)); //allocate array m
	for (i = 0; i < dim_y; i++) {
		(*m)[i] = malloc(dim_x * sizeof(double));
	}

}
