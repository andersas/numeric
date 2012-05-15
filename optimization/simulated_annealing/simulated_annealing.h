#ifndef SIMANNEAL
#define SIMANNEAL
#include <stdio.h>

double simulated_annealing(size_t dimension, double *result, \
			          double (*f)(const double x[dimension]), \
					const double lower_limits[], \
					const double higher_limits[], \
					double max_step_distance, \
					int n_runs, double gamma, int n_steps);

#endif
