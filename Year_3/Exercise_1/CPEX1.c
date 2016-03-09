/*
 ============================================================================
 Name        : CPEX.c
 Author      : Sam Low

 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double DET(double **a, int n); //function prototype
void cofactor(double **a, int n, double ***m); //a wild triple pointer has appeared
void transpose(double **a, int n);

/*Program will request dimension n, then print an nxn matrix of
 random numbers to screen, followed by its inverse */

int main() {
	int n, j, i; //n is matrix dimension, i,j just counters
	double **A; //pointer to pointer to double
	double **cof;

	printf("Enter the number of rows and columns n:");
	scanf("%d", &n);

	A = malloc(n * sizeof(double*)); //allocating memory to 1D array

	for (i = 0; i < n; i++) {
		A[i] = malloc(n * sizeof(double)); //allocate to each element, spanning into 2D
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			A[i][j] = (20 * rand()) / (double) RAND_MAX;
		}
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("A[%d][%d] = %lf\n", i, j, A[i][j]);
		}
	}
	printf("\n");

	cofactor(A, n, &cof);
	/* finds the cofactor matrix, and stores physically over the location of cof
	 in memory. should be able to read out values 'from cof' */

	transpose(cof, n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			cof[i][j] = (1 / DET(A, n)) * cof[i][j]; //cof is now the inverse of A
		}
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("Inv[%d][%d] = %lf\n", i, j, cof[i][j]);
		}
	}

}

/*
 Determinant of A:
 1. Run across top of A. i.e A00, A01, A02... etc
 2. At each element A_0j copy all elements NOT in rows/cols i=0,j into a new matrix.
 3. Find the determinant of this 'minor' matrix
 4. Calculate DetA using expansion of minors. ie SUM: (-1)^(i+j) A[i][j]*Det(M_ij)
 */

double DET(double **a, int n) {
	int i, col_m, col_A, col_m2;
	double det = 0;
	double **m = NULL;	//new matrix. will be manipulated as each minor in turn.

	/*if the matrix is 2x2, calculate the determinant. if it is bigger want to use expansion of minors
	 to reduce it over and over again into a number of smaller matrices: minors. eventually reach 2x2
	 and calculate determinant. do this as many times as necessary.*/

	if (n == 2) { 		//if the matrix reduces to 2x2, calc det and finish.
		det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	}

	else {

		/*
		 col_A = SOURCE (A) COLUMN COUNTER
		 col_m = MINOR COLUMN COUNTER
		 col_m2 = MINOR COLUMN COUNTER (COPY TARGET)
		 */

		det = 0;

		for (col_A = 0; col_A < n; col_A++) //working across top of A
				{

			m = (double **) malloc((n - 1) * sizeof(double *)); //memory for n-1 *doubles

			for (i = 0; i < n - 1; i++)
				m[i] = (double *) malloc((n - 1) * sizeof(double)); //memory for doubles at each *double

			/*now BUILD THE MINOR MATRIX excluding col_A and every row in turn.
			 reduce what remains into a space 1 row/col smaller.*/

			for (i = 1; i < n; i++) {
				col_m2 = 0; // copy forward the source matrix, excluding one column

				for (col_m = 0; col_m < n; col_m++) {
					if (col_m != col_A) { // if col_m matches col_A DON'T copy into minor
						m[i - 1][col_m2] = a[i][col_m]; // otherwise copy src element into minor [i-1][col_m2]

						col_m2++;
					}
				}
			}

			det += pow(-1.0, col_A + 2.0) * a[0][col_A] * DET(m, n - 1); //EXPANSION OF MINORS FORMULA

			for (i = 0; i < n - 1; i++) {
				free(m[i]);
			}	// free the storage allocated to this minor's set of pointers
			free(m);     // free the storage for the original pointer to pointer

		}
	}
	return (det);
}

/*
 Cofactor matrix of A:
 1. run through matrix A
 2. at each element A_ij copy all elements NOT in rows/cols i,j into a new matrix.
 3. Find the determinant of this matrix
 4. Copy into cofactor matrix
 5. Tranpose it
 */

void cofactor(double **a, int n, double ***m) {
	int Ai = 0, Aj = 0;
	int i = 0, j = 0, Ci = 0, Cj = 0;
	double **c, factor;

//create cofactor matrix cof to store the determinants of matrices 'c'
	(*m) = (double **) malloc((n) * sizeof(double *));
	for (i = 0; i < n; i++)
		(*m)[i] = (double *) malloc((n) * sizeof(double));

//scanning through every element of A
	for (Ai = 0; Ai < n; Ai++) {
		for (Aj = 0; Aj < n; Aj++) {

			//at each elemment A[Ai][Aj] create minor matrix c
			c = (double **) malloc((n - 1) * sizeof(double *));
			for (i = 0; i < n - 1; i++)
				c[i] = (double *) malloc((n - 1) * sizeof(double));

			//at each element A[Ai][Aj], scan through elements A[i][j] again
			Ci = 0;
			for (i = 0; i < n; i++) {
				Cj = 0;
				if (i != Ai) {
					for (j = 0; j < n; j++) {

						//if element is NOT in same row/col as A[Ai][Aj]
						if (j != Aj) {

							c[Ci][Cj] = a[i][j]; //copy it to matrix c
							Cj++;
						}
					}
					Ci++;
				}
			}
			//matrix c should be complete now
			factor = pow(-1, Ai + Aj + 2) * DET(c, n - 1);
			(*m)[Ai][Aj] = factor; 	//matrix m now holds the cofactors.
			//This does NOT mean that matrix cof in main does.

		}
	}		//finish looping through Ai and Aj

}		//f.end

void transpose(double **a, int n) {
	int i, j;
	double x;
	/*need x as a temporary store for cof_ij. Just writing a[i][j]=a[j][i] in the loop
	 is insufficient. you cannot SWAP using this expression, only change the LHS. to SWAP
	 you need two distinct expressions: put A in B and put what WAS in B in A*/

	for (i = 1; i < n; i++) {
		for (j = 0; j < i; j++) {

			x = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = x;
		}
	}

}
