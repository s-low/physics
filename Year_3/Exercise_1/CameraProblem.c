/*
 ============================================================================
 Name        : CameraProblem.c
 Author      : Sam Low
 ============================================================================
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

/*Program will request cable identity, then print to file:
 displacement from cable's attachment point, tension in cable */

int main(void) {
	int s,c;
	double x1, y1, x2, y2, x3, y3, mod1, mod2, mod3;
	double base = 121.24355653; //base of the triangle
	FILE *datafile;
	double dx; //increment in space x,y

	printf("Which cable: 1, 2, or 3?");
	scanf("%d", &c);
	c=c-1;

	/*
	 |x1 x2 x3||T1|		| 0|
	 |y1 y2 y3||T2|	=	| 0|
	 |z1 z2 z3||T3|		|mg|

	 A(mat) T(vec) = B(vec)
	 */

	gsl_matrix * A = gsl_matrix_alloc(3, 3);
	gsl_vector * T = gsl_vector_alloc(3);
	gsl_vector * F = gsl_vector_alloc(3);

	gsl_permutation * p = gsl_permutation_alloc(3); //permutation matrix. required for LU decomp.

	dx = 0.5; //increment of 50cm, sensible for most plotting software but 0.1 is better

	datafile = fopen("data.txt", "w");
	if (datafile != NULL) {
		for (x1 = 0; x1 < base; x1 += dx) { //scans square containing triangle
			for (y1 = 0; y1 < base; y1 += dx) {

				//now omit points inside the square which do not lie inside the triangle
				//i.e bounded by two straight lines at edges of triangle
				if ((y1 < x1 * sqrt(3)) && (y1 < (210 - x1 * sqrt(3)))) {

					x2 = x1 - (base); //calc the other displacement coords rel. to corners 2,3
					y2 = y1;

					x3 = x1 - (base / 2);
					y3 = (-105) + y1;

					mod1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(10, 2));//modulus of each the displacement
					mod2 = sqrt(pow(x2, 2) + pow(y2, 2) + pow(10, 2));
					mod3 = sqrt(pow(x3, 2) + pow(y3, 2) + pow(10, 2));

					gsl_matrix_set(A, 0, 0, x1 / mod1);	//fill A with unit vector components
					gsl_matrix_set(A, 0, 1, x2 / mod2);
					gsl_matrix_set(A, 0, 2, x3 / mod3);
					gsl_matrix_set(A, 1, 0, y1 / mod1);
					gsl_matrix_set(A, 1, 1, y2 / mod2);
					gsl_matrix_set(A, 1, 2, y3 / mod3);
					gsl_matrix_set(A, 2, 0, 10 / mod1);
					gsl_matrix_set(A, 2, 1, 10 / mod2);
					gsl_matrix_set(A, 2, 2, 10 / mod3);

					gsl_vector_set(F, 0, 0);		//fill F with Forces 0,0,mg
					gsl_vector_set(F, 1, 0);
					gsl_vector_set(F, 2, (50 * 9.81));

					gsl_linalg_LU_decomp(A, p, &s);			//LU decomp PA = LU
					gsl_linalg_LU_solve(A, p, F, T);	//A is U at this point

					fprintf(datafile, "%lf\t%lf\t%lf\t\n", x1, y1,
							gsl_vector_get(T, c)); // prints to file: X,Y,Tension in cable1
				}
			}
		}
	}

	fclose(datafile);
}

