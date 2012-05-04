#include <stdlib.h>
#include <string.h>
#include "data.h"
#include "qspline.h"

double spline_evaluate(struct qsplines_workspace *state, double z);

// Calculate the coefficients a_i and store them in
// a qsplines_workspace together with a method to evaluate
// the relevant spline at a point.

struct qsplines_workspace * qspline_init(int npoints, const struct point points[]) {

	struct qsplines_workspace *state;
	double *a, *b;
	double *dy,*dx;
	int i;

	state = (struct qsplines_workspace *) malloc(sizeof(struct qsplines_workspace));

	a = (double *) malloc((npoints-1)*sizeof(double));
	b = (double *) malloc((npoints-1)*sizeof(double));
	dy = (double *) malloc((npoints-1)*sizeof(double));
	dx = (double *) malloc((npoints-1)*sizeof(double));


	// Calculate differences between points
	
	for (i = 0; i < npoints -1; i++) {
		dy[i] = points[i+1].y - points[i].y;
		dx[i] = points[i+1].x - points[i].x;
	}

	// Set last and first a to zero and run both recursions
	// then average the result
	
	a[0] = 0;
	b[npoints-2] = 0;

	for (i = 0; i < npoints - 2; i++) {
	
		a[i + 1] = (dy[i+1]/dx[i+1] - dy[i]/dx[i] - a[i]*dx[i])/dx[i+1];

		b[npoints - i - 3] = (dy[npoints - i - 2]/dx[npoints - i - 2] \
		- dy[npoints - i - 3] / dx[npoints -i - 3] \
		- b[npoints - i - 2] * dx[npoints - i - 2]) / dx[i];

	}

	// Average
	
	for (i = 0; i < npoints - 1; i++) {
		a[i] = (a[i] + b[i])/2.0;
	}
	

	free(b); // Don't need this anymore
	free(dy);
	free(dx);
	
	// Generate return structure

	state->a = a;
	state->count = npoints;
	state->eval = &spline_evaluate;
	state->last_index = 0;
	state->data = (struct point *) malloc(sizeof(struct point) * npoints);

	memcpy((void *) state->data, (void *) points, sizeof(struct point) * npoints);

	return state;
}

// This reduces typing:
#define x(i) state->data[i].x
#define y(i) state->data[i].y

double spline_evaluate(struct qsplines_workspace *state, double z) {

	int min, max, try, i;
	int last_result;
	int npoints;

	last_result = state->last_index;
	npoints = state->count;

	// See if we can be clever and guess what spline
	// we need based on what spline we used last.
	// We can avoid doing a binary search each call with this
	// if the user repeatedly asks for a value slightly
	// larger than the previous one like he would if
	// he was plotting the function.
	
	if (x(last_result) <= z && z <= x(last_result + 1)) {
		min = last_result;
	} else {
		if (last_result < npoints - 1 && x(last_result + 1) <= z && z <= x(last_result + 2)) {
                        min = state->last_index += 1;
                } else if (last_result > 0 && (x(last_result - 1) <= z && z <= x(last_result))) {
                        min = state->last_index -= 1;
                } else {

			// Clever guessing failed, resort to
			// binary search
		
			min = 0;
			max = npoints - 1;
			while (min + 1 < max) {
				try = (max + min)/2;
				if (x(try) > z)
					max = try;
				else
					min = try;
			}
			state->last_index = min;
		}
	}

	// Now evaluate spline i in z and get out of here
	
	i = min;

	return y(i) + ((y(i+1) - y(i)) / (x(i+1) - x(i))) * (z-x(i)) + \
	       state->a[i]*(z-x(i))*(z-x(i+1));

}


