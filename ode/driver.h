#ifndef DRIVER
#define DRIVER


typedef struct driver_r {
	double x;
	double *y;
	struct driver_r *next;
} driver_r;


void driver_r_destroy(driver_r *);


driver_r * driver(void (*f)(double x, const double y[], double result[]), int dim, double start, double stop, double y0[dim], double acc, double eps, double h, int method);

#endif

