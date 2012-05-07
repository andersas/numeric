#include <stdio.h>
#include "../linalg/linalg.h"
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include <math.h>


matrix * get_initial_point(size_t dim, matrix *lower, matrix *higher, gsl_rng *rng) {

	
	


}


void annealing_step(size_t dim matrix *x, double step_size, gsl_rng *rng) {
// calculate a random direction and add it to x
	
	int i;
	double len;
	matrix * dx = matrix_create(dim,1);

	do {
		len = 0.0;	
		for (i = 0; i < dim; i++) {
			dx->a[i][0] = gsl_rng_uniform_pos(rng);
			len += (dx->a[i][0])*(dx->a[i][0]);
		}
		len = sqrt(len);
	} while (len > 1.0); // Get something not outside the unit sphere

	for (i = 0; i < dim; i++) {
		dx->a[i][0] /= len;
		dx->a[i][0] *= step_size*gsl_rng_uniform(rng);
		x->a[i][0] += dx->a[i][0];
	}

}


matrix * simulated_annealing(size_t dimension, double (*f)(const double x[dimension]), const matrix *lower_limits, const matrix *higher_limits, double max_step_distance, int n_runs, double temperature, double dT) {

	
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



	}


	fclose(handle);
}

