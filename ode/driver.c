#include "step.h"
#include "driver.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>


// Returns driver_r, a linked list of points,
// and the last point in stop, y0.
driver_r * driver(void (*f)(double x, const double y[], double result[]), int dim, double start, double stop, double y0[dim], double acc, double eps, double h, int method) {

	driver_r *res, *last, *next;
	double y[dim];
	double x = start;
	double tol;
	double err;
	double Power = 0.25;
	double Safety = 0.95;
	int i;

	double (*step)(void (*f)(double x, const double y[], double result[]), int dim, double x0, double y0[dim], double h);

	if (method == 1)
		step = &rkstep1;
	else
		step = &rkstep4;

	memcpy(&y, y0, dim*sizeof(double));


	// Initialize first node in linked list.
	res = (driver_r *) malloc(sizeof(driver_r));
	last = res;
	last->x = x;
	last->y = (double *) malloc(sizeof(double)*dim);
	memcpy(last->y, y, sizeof(double)*dim);
	last->next = NULL;

	while (x < stop) {

	err = step(f, dim, x, y, h);
	tol = 0.0;
	for (i = 0; i<dim; i++) {
		tol += y[i]*y[i];
	}
	tol = (eps*sqrt(tol) + acc)*sqrt(h/(stop-start));

	if (tol > err) {
		// Accept the point, advance to next iteration
		x += h;
		next = (driver_r *) malloc(sizeof(driver_r));
		next->x = x;
		next->y = (double *) malloc(dim*sizeof(double));
		memcpy(next->y, y, dim*sizeof(double));
		next->next = NULL;
		last->next = next;
		last = next;
	} else {
		memcpy(y,last->y,dim*sizeof(double)); // Restore previous y values
	}
	if (err > tol/100) {
		h = h * pow(tol/err,Power) * Safety;
	} else {
		h *= 1.05; // If error is too small, just add
			   // 5 % to next step size.
	}
	
	}

	memcpy(y0, y, dim*sizeof(double)); // Return y through y0

	return res;

}


void driver_r_destroy(driver_r *r) {

	driver_r *next;

	while (r != NULL) {
		next = r->next;
		free(r->y);
		free(r);
		r = next;
	}

}

