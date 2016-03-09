/*
 ============================================================================
 Name        : RejectAccept.c
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
	double y, x, chi;
	float binsize;
	int j, bin, n, k;
	float *m;

	printf("Range for random numbers is ~3.2\n\nEnter the number of bins:");
	scanf("%d", &n);

	printf("\nk = 10^");
	scanf("%d", &k);

	k = pow(10, k);

	binsize = 3.2 / n;

	m = malloc(n * sizeof(float)); //1D array to store the count for each bin

	datafile = fopen("RNG_RejAcc.txt", "w");
	if (datafile != NULL) //check to see if the file was actually opened
	{
		gsl_rng * r;
		r = gsl_rng_alloc(gsl_rng_taus);

		for (j = 0; j < k; j++) {

			double y = gsl_rng_uniform(r);
			double x = gsl_rng_uniform(r);

			x = M_PI * x; //increase range to 0->pi

			if (y < sin(x)) { //accept if less than *test(y)*

				x = (x / binsize); //reduce x to multiple of bins

				m[(int) x] = m[(int) x] + 1;
			}

		}
		gsl_rng_free(r);

		//FILE FORMAT: x, frequency
		//note frequencies are normalised to fit under a HALF HEIGHT sin curve

		for (bin = 0; bin < n; bin += 1) {
			fprintf(datafile, "%.3f\t%f\n", (bin * binsize) + (binsize / 2),
					(float) (n * m[bin]) / (2 * k));
		}

		chi = chisquare(n, m, k, binsize);
		printf("chi square = %lf", chi); //puts chisquare stat to screen

		printf("\n\nfile written\n");

	}
	fclose(datafile);
	return 0;
}

double chisquare(int n, float *m, int k, float binsize) {
//subroutine calculates chisquare test for fitting

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
	} //normalise all data to range 0->1 (just for purposes of test)

	for (bin = 0; bin < n; bin += 1) {
		mean = mean + m[bin];
	} //calculate mean bin value (sum all values)

	mean = mean / n; //finish mean calculation

	for (bin = 0; bin < n; bin += 1) {
		var_i = m[bin] - mean;
		var = var + pow(var_i, 2);
	} //calculate variance

	var = var / n; //finish calculation

	for (bin = 0; bin < n; bin += 1) {
		chi_i = pow((sin(bin * binsize) - m[bin]), 2);
		chi_i = chi_i / var;
		chi = chi + chi_i;
	} //calculate total chi^2

	return chi;

}

