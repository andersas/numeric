#include <stdio.h>
#include "leastSquares.h"
#include "../linalg/linalg.h"

matrix *fittocircle(char *filename, double *x0, double *y0, double *radius);

int main (int argc, char **argv) {

	double x0, y0, radius;
	matrix *cov;

	// Read circ.bat, generated by gencircle.py.
	// This is a circle with radius 15.2 with center 2.4, -3.6
	// and radius 15.2 with gaussian noise.

	cov = fittocircle("circ.dat", &x0, &y0, &radius);

	printf("Circle has center (%f, %f) with radius %f.\n", x0, y0, radius);

	// See circfit.png

	matrix_store(cov, "cov.out");

	// There are great covariances in the third parameter,
	// as this is a direct function of the two other parameters.

	matrix_free(cov);

	


return 0;
}