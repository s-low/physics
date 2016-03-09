/*
 ============================================================================
 Name        : GaussSeidel.c
 Author      : Sam Low
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*edges is used to set all edges of matrix to same value
 * allocate is used to allocate a 2d array*/

void edges(int dim_x, int dim_y, double ***m, double value);
void allocate(int dim_x, int dim_y, double ***m);

int main() {

	int i, i1 = 0, j, loop;
	double **V, **Vstore; //jacobi method - two complete copies of grid
	int dimx, dimy; //grid dimensions
	double diff, perc, X;
	FILE *datafile;

	printf("Program uses plate capacitor with extremal corner of y=60cm, x=60cm\n");
	printf("Be sure to exceed these dimensions\n");

	printf("\nx dimension in whole centimetres:");
	scanf("%d", &dimx);

	printf("\ny dimension in whole centimetres:");
	scanf("%d", &dimy);

	printf("\nEnter the the convergence criterion in percent:");
	scanf("%lf", &X);
	//convergence criterion. Stop iterating if no node change by more than X%

	allocate(dimx, dimy, &V);
	allocate(dimx, dimy, &Vstore);

	for (i = 0; i < dimy; i++) { //fill V with initial guess
		for (j = 0; j < dimx; j++) {
			V[i][j] = 0;
		}
	}

	/*BOUNDARY CONDITIONS*/

	for (i = 40; i < 60; i++) { //+ve plate
		V[i][60] = 10;
	}

	for (i = 40; i < 60; i++) { //-ve plate
		V[i][40] = -10;
	}

	edges(dimx, dimy, &V, 0);

	/*Boundaries set, now iterate finite difference method*/

	datafile = fopen("GaussSeidel.txt", "w");
	if (datafile != NULL) //check to see if the file was actually opened
	{

		loop = 1;
		while (loop == 1) { //implementing the convergence criterion as a while loop

			loop = 0; //unless this is reset to 1 below, the method won't repeat
			for (i = 1; i < dimy - 1; i++)
				for (j = 1; j < dimx - 1; j++) {

					Vstore[i][j] = V[i][j]; //save old V values somewhere

					V[i][j] = 0.25
							* ((V[i + 1][j]) + (V[i - 1][j]) + (V[i][j + 1])
									+ (V[i][j - 1]));
				}

			/*the key difference between GS and J method is above:
			 new V is found from within V itself */

			for (i = 40; i < 60; i++) { //+ve plate
				Vstore[i][60] = 10;
				V[i][60] = 10;
			}

			for (i = 40; i < 60; i++) { //-ve plate
				Vstore[i][40] = -10;
				V[i][40] = -10;
			}

			edges(dimx, dimy, &V, 0); //reset all edges to zero
			edges(dimx, dimy, &Vstore, 0);

			for (i = 1; i < dimy - 1; i++) //check every element for change
				for (j = 1; j < dimx - 1; j++) {

					diff = pow(V[i][j] - Vstore[i][j], 2); //magnitude
					perc = 100 * (sqrt(diff) / ((V[i][j] + Vstore[i][j]) / 2));

					if (perc > X) { //convergence criterion: no node changes > X%
						loop = 1;
					}
				}

		}

		for (i = 0; i < dimy; i++)
			for (j = 0; j < dimx; j++) {

				if (i != i1) {
					fprintf(datafile, "\n"); //quirk of using gnuplot with pm3d
				}
				i1 = i;

				fprintf(datafile, "%d\t%d\t%lf\n", i, j, V[i][j]);

			}
	}

	fclose(datafile);
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

