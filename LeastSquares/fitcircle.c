#include <stdio.h>
#include <math.h>
#include "leastSquares.h"
#include "../linalg/linalg.h"


/* Fit a number of data points (xi, yi) to a circle
 * (x-x_0)^2 + (y-y_0)^2 = r^2,
 * or
 * x^2 + y^2 = (r^2 - x_0^2 - y_0^2) + 2xx_0 + 2yy_0.
 *
 * If we call z_0 = (r^2 - x_0^2 - y_0^2) we get
 * 2xx_0 + 2yy_0 + z_0 = x^2 + y^2, and thus,
 * f1(i) = 2xi
 * f2(i) = 2yi
 * f3(i) = 1,
 * we need to find x_0, y_0, z_0 which minimizes the error in
 * x_0f1(i) + y_0f2(i) + z_0f3(i) = x^2 + y^2 */

matrix *xs, *ys, *dys;

double f1(double x) {
	return 2*xs->a[lrint(x)][0];
}

double f2(double x) {
	return 2*ys->a[lrint(x)][0];
}

double f3(double x) {
	return 1.0;
}

matrix *fittocircle(char *filename, double *x0, double *y0, double *radius) {

	matrix *data;
	matrix *x, *y;
        matrix *coeff, *cov;
	double (*f[])(double) = {&f1, &f2, &f3};
	size_t i,n,m;

	data = matrix_load(filename);
	n = data->n;
	m = data->m;
	if (m != 3) {
		fprintf(stderr,"Invalid amount of data columns.\n");
		exit(-1);
	}
	
	xs = matrix_create(n,1);
	x = matrix_create(n,1);
	ys = matrix_create(n,1);
	y = matrix_create(n,1);
	dys = matrix_create(n,1);

	for (i = 0; i < n; i++) {
		xs->a[i][0] = data->a[i][0];
		ys->a[i][0] = data->a[i][1];
		dys->a[i][0] = data->a[i][0];
		x->a[i][0] = i;
		y->a[i][0] = xs->a[i][0]*xs->a[i][0] + ys->a[i][0]*ys->a[i][0];
	}
	matrix_free(data);

	coeff = matrix_create(3,1);
	cov = matrix_create(3,3);

	leastSquares(n,3,x,y,dys,f,coeff,cov);

	*x0 = coeff->a[0][0];
	*y0 = coeff->a[1][0];
	*radius = sqrt(coeff->a[2][0] + coeff->a[0][0]*coeff->a[0][0] + coeff->a[1][0]*coeff->a[1][0]);

	free(xs);
	free(ys);
	free(dys);
	free(x);
	free(y);
	free(coeff);

	return cov;	
}
