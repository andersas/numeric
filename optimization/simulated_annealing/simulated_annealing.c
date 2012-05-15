#include <stdio.h>
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
			dx[i] = 2*(gsl_rng_uniform_pos(rng) - 0.5);
			len += (dx[i])*(dx[i]);
		}
		len = sqrt(len);
	} while (len > 1.0); // Get something not outside the unit sphere

	for (i = 0; i < dim; i++) {
		dx[i] /= len;
		dx[i] *= step_size*gsl_rng_uniform(rng);
		x[i] += dx[i];
	}

printf("dx: %f\n", dx[0]);
}

__inline__ void accept_step(size_t dimension, double *bestf, double bestx[dimension], double *fnow, double *fnext, double x[dimension],double candidate[dimension]) {
			memcpy(x,candidate,sizeof(double)*dimension);
			*fnow = *fnext;
			if (*bestf > *fnext) {
				memcpy(bestx,x,sizeof(double)*dimension);
				*bestf = *fnext;
			}

			printf("Accepted new step: x: %f, f(x) = %f.\n", x[0],*fnow);
}

void ignore_step(size_t dimension, double new[dimension], double old[dimension]){
	memcpy(new,old,sizeof(double)*dimension);
}

double simulated_annealing(size_t dimension, double *result, \
			          double (*f)(const double x[dimension]), \
					const double lower_limits[], \
					const double higher_limits[], \
					double max_step_distance, \
					int n_runs, double gamma, int n_steps) {

	
	double T;
	double x[dimension];
	double z[dimension];
	double bestx[dimension];
	double bestf;
	double fx,fz;
	double probability;
	unsigned long int seed;
	int step;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	FILE *handle = fopen("/dev/urandom", "r");

	fread(&seed, sizeof(unsigned long int), 1, handle);
	gsl_rng_set(rng,seed);


	get_initial_point(dimension,x,lower_limits,higher_limits,rng);
	
	memcpy(z,x,sizeof(double)*dimension);
	memcpy(bestx,x,sizeof(double)*dimension);
	bestf = f(x);
	fx = bestf;

	step = 0;
	do {
		
		annealing_step(dimension, z, max_step_distance, rng);

		fz = f(z);

		if (fz <= fx) {
			accept_step(dimension,&bestf,bestx,&fx,&fz,x,z);
		} else {
			ignore_step(dimension,z,x);
		/*	T = gamma/log(step+2);
			probability = exp((fz-fx)/T);
			if (gsl_rng_uniform_pos(rng) < probability) {
				accept_step(dimension,&bestf,bestx,&fx,&fz,x,z);
			} else {
				ignore_step(dimension,z,x);
			}*/
		}
		
		step++;	
	} while (step < n_steps);

	
	fclose(handle);

	memcpy(result,bestx,sizeof(double)*dimension);
	return bestf;
}

