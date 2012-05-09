#include <stdio.h>
#include "../linalg/linalg.h"
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include <string.h>
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


double simulated_annealing(size_t dimension, double *result, \
			          double (*f)(const double x[dimension]), \ 
					const double lower_limits[], \
					const double higher_limits[], \
					double max_step_distance, \
					int n_runs, double temperature) {

	
	double T;
	double x[dimension];
	double z[dimension];
	double bestx[dimension];
	double bestf;
	double fx,fz;
	unsigned long int seed;
	int step;
	gsl_rng *rnd = gsl_rng_alloc(gsl_rng_mt19937);
	FILE *handle = fopen("/dev/urandom", "r");

	fread(&seed, sizeof(unsigned long int), 1, handle);
	gsl_rng_set(rnd,seed);


	get_initial_point(dimension,x,lower_limits,higher_limits,rng);
	
	memcpy(z,x,sizeof(double)*dimension);
	memcpy(bestx,x,sizeof(double)*dimension);
	bestf = f(x);
	fx = bestf;

	step = 0;
	do {
		
		annealing_step(dimension, z, max_step_size, rng);

		fz = f(z);

		if (fz <= fx) {
			
		} else {

		}
		
		
		//   t = gamma/log(k+2), k = number of steps,
		//	gamme problem constant

		//	probabilities = min(1,exp(f(next)-f(now)/T));

	} while();


	fclose(handle);
}

