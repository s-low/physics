/*
 ============================================================================
 Name        : GammaRay.c
 Author      : Sam Low
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <gsl/gsl_randist.h>

void allocate(int dim_x, int dim_y, double ***m);
//allocates memory to a 2D array of doubles x by y

int main(void) {

	double mu, x, y, z, rad, phi, theta, binwidth, size, sigmax, sigmay;
	int i, j, k, n, binx, biny;
	double **m;
	FILE *datafile;

	printf("Detector array is 2m x 2m\n\n");
	printf("n = number of pixels (i.e nxn)\n");
	printf("\nn = ");
	scanf("%d", &n); //array pixels nxn

	size = 2; //2m by 2m array

	//standard dev = FWHM(resolution) / 2.3548)
	sigmay = 0.3 / 2.3548;
	sigmax = 0.1 / 2.3548;

	printf("k = total random numbers to generate\n");
	printf("\nk = 10^");
	scanf("%d", &k);

	k = pow(10, k); //number of random numbers to generate
	mu = 1.04; //mean path length = 1.04m
	binwidth = size / n; //width of a pixel in m
	allocate(n, n, &m); //nxn pixel array

	datafile = fopen("gamma_simulation.txt", "w");
	if (datafile != NULL) {

		gsl_rng * r;
		r = gsl_rng_alloc(gsl_rng_taus);

		for (j = 0; j < k; j++) {

			gsl_ran_dir_3d(r, &x, &y, &z); //generate 3d vector of MAG:1

			//use spherical coords to get direction in angles
			//take care with phi, cannot always just take arctan of y/x

			theta = acos(z);

			if (x > 0) {
				phi = atan(y / x);
			}

			else if (x == 0) {
				phi = M_PI * (y / fabs(y));
			}

			else if (x < 0) {
				phi = atan(y / x) + (M_PI * (y / fabs(y)));
			}

			double d = gsl_ran_exponential(r, mu);
			//exponential distributed random path length

			rad = (2 - d) / cos(theta); //extrapolate direction to screen
			x = rad * sin(theta) * cos(phi); //x and y coords on screen
			y = rad * sin(theta) * sin(phi);

			x = x + gsl_ran_gaussian(r, sigmax); //add gaussian smearing from resolution
			y = y + gsl_ran_gaussian(r, sigmay);

			//if statement to check if the coordinate in on the detector screen
			if ((fabs(x) < size) && (fabs(y) < size)) {

				binx = (int) ((x * (n / size)) + (n / 2));
				biny = (int) ((y * (n / size)) + (n / 2));
				// n/2 moves origin to centre of screen

				if ((binx < n) && (biny < n) && (binx > 0) && (biny > 0)) {
					m[binx][biny]++;
				}
				//increment the frequency at the bin (if the bin exists)

			}

		}
		gsl_rng_free(r);

		//FILE FORMAT: x coord, y coord, no. of events
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				fprintf(datafile, "%lf\t%lf\t%lf\n", i * binwidth, j * binwidth,
						(m[i][j]));
			}
			fprintf(datafile, "\n");
		}

		printf("\n\nfile written\n");
	}
	fclose(datafile);
}

void allocate(int dim_x, int dim_y, double ***m) { //function to allocate 2d array of doubles
	int i; //just a counter

	(*m) = malloc(dim_y * sizeof(double*)); //allocate array m
	for (i = 0; i < dim_y; i++) {
		(*m)[i] = malloc(dim_x * sizeof(double));
	}

}

