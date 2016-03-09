/*
 ============================================================================
 Name        : RandomNumbers1.c
 Author      : Sam Low
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>

double chisquare(int n, float *m, int k, float binsize);
//subroutine calculates the value of chisquared fitted against sin(x)

int main(void) {
	FILE *datafile;
	double xreq, chi;
	int j, bin, n, k;
	float *m, binsize;

	printf("Range for random numbers is ~3.2\n\nEnter the number of bins:");
	scanf("%d", &n);

	printf("\nk = 10^");
	scanf("%d", &k);

	k = pow(10, k);

	binsize = 3.2 / n;

	m = malloc(n * sizeof(float)); //use 1D array to store the count for each bin

	datafile = fopen("RNG_Analytic.txt", "w");
	if (datafile != NULL) //check to see if the file was actually opened
	{
		gsl_rng * r;
		r = gsl_rng_alloc(gsl_rng_taus);

		for (j = 0; j < k; j++) {

			double xgen = gsl_rng_uniform(r); //generate uniform random x
			xreq = acos((1 - 2 * xgen)); //make transform (see report)

			xreq = (xreq / binsize);
			//reduces random number to a multiple of binwidths

			m[(int) xreq] = m[(int) xreq] + 1;
			//cast mutiple of binwidths to int, and increment corresp. bin
		}
		gsl_rng_free(r);


		//FILE FORMAT: x, frequency
		//note frequencies are normalised to fit under a HALF HEIGHT sin curve
		for (bin = 0; bin < n; bin += 1) {
			fprintf(datafile, "%.3f\t%f\n", (bin * binsize) + (binsize / 2),
					(float) (n * m[bin]) / (M_PI * k));
		}

		chi = chisquare(n, m, k, binsize);
		printf("chi square = %lf", chi); //puts chisquare stat to screen

		printf("\n\nfile written\n");

	}
	fclose(datafile);
}

double chisquare(int n, float *m, int k, float binsize) {
//calculates chisquared test

	double var = 0, var_i = 0, chi = 0, chi_i = 0;
	double mean = 0;
	int bin, max = 0;

	for (bin = 0; bin < n; bin += 1) {
		if (m[bin] > max) {
			max = m[bin];
		}
	} //find max frequency in range by scanning through and replacing if bigger

	for (bin = 0; bin < n; bin += 1) {
		m[bin] = m[bin] / max;
	} //normalise all data to range 0->1 (only for purposes of test)

	for (bin = 0; bin < n; bin += 1) {
		mean = mean + m[bin];
	} //calculate mean bin value

	mean = mean / n; //finish calculation of mean

	for (bin = 0; bin < n; bin += 1) {
		var_i = m[bin] - mean;
		var = var + pow(var_i, 2);
	} //calculate variance

	var = var / n; //finish calc.

	for (bin = 0; bin < n; bin += 1) {
		chi_i = pow((sin(bin * binsize) - m[bin]), 2);
		chi_i = chi_i / var;
		chi = chi + chi_i;
	} //calculate total chi^2

	return chi;

}

