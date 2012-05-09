#include <stdio.h>
#include "../linalg/linalg.h"
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include <math.h>


// Writes initial point to *initial
void get_initial_point(size_t dim, double *initial, const double lower[], const double higher[], gsl_rng *rng) {

	
	int i;

	for (i = 0; i < dim; i++) {
		initial[i] = gsl_rng_uniform_pos(rng)*(higher[i]-lower[i]) - lower[i];
	}
}


void annealing_step(size_t dim, double *x, double step_size, gsl_rng *rng) {
// calculate a random direction and add it to x
	
	int i;
	double len;
	double dx[dim];

	do {
		len = 0.0;	
		for (i = 0; i < dim; i++) {
			dx[i][0] = gsl_rng_uniform_pos(rng);
			len += (dx[i][0])*(dx[i][0]);
		}
		len = sqrt(len);
	} while (len > 1.0); // Get something not outside the unit sphere

	for (i = 0; i < dim; i++) {
		dx[i][0] /= len;
		dx[i][0] *= step_size*gsl_rng_uniform(rng);
		x[i][0] += dx[i][0];
	}

}


void simulated_annealing(size_t dimension, double *result, \
			          double (*f)(const double x[dimension]), \ 
					const double lower_limits[], \
					const double higher_limits[], \
					double max_step_distance, \
					int n_runs, double temperature, double dT) {

	
	double T;
	double bestx[dimension];
	double bestf;
	unsigned long int seed;
	gsl_rng *rnd = gsl_rng_alloc(gsl_rng_mt19937);
	FILE *handle = fopen("/dev/urandom", "r");

	fread(&seed, sizeof(unsigned long int), 1, handle);
	gsl_rng_set(rnd,seed);

	T = temperature;
	while (T > 0) {

		
		   t = gamma/log(k+2), k = number of steps,
			gamme problem constant

			probabilities = min(1,exp(f(next)-f(now)/T));

	}


	fclose(handle);
}

