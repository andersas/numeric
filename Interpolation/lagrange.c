#include "lagrange.h"


/* Evaluate the lagrange polynomial for npoints points pointed
 * to by *points at a point z.
 */

double lagrange_interpolate(int npoints, struct point *points, double z) {

	int i, k;
	double numerator, denominator, sum;

	/* Calculate the interpolated value directly from the lagrange form */

	sum = 0.0;
	for (i = 0; i < npoints; i++) {

		numerator = denominator = 1.0;
		for (k = 0; k < npoints; k++) {
			if (k != i) {
				numerator *= z - points[k].x;
				denominator *= points[i].x - points[k].x;
			}
		}
		sum += points[i].y * numerator/denominator;
	}

	return sum;

}



