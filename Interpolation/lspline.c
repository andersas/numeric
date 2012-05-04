#include <stdlib.h>
#include "data.h"
#include "lspline.h"

// This saves a lot of typing. Just remember to use () instead of []
#define x(i) points[i].x
#define y(i) points[i].y

// NOT THREAD SAFE! (we use a static variable)
double lspline_interpolate(int npoints, struct point *points, double z) {

	int min,max, try, i;
	static int last_result = 0;  // skip the binary search
				     // when evaluating
				     // along the same spline as last time

	if (x(last_result) <= z && z <= x(last_result + 1)) {
		min = last_result;
	} else {
		if (last_result < npoints - 1 && x(last_result + 1) <= z && z <= x(last_result + 2)) {
			min = last_result += 1;
		} else if (last_result > 0 && (x(last_result - 1) <= z && z <= x(last_result))) {
			min = last_result -= 1;
		} else {
			min = 0;
			max = npoints - 1;
			// Figure out what spline we need to use
			// by binary search for i such that x[i] <= z <= x[i+1] 
			while (min + 1 < max) {
				try = (max+min)/2;
				if (x(try) > z)
					max = try;
				else
					min = try;
			}
			last_result = min;
		}
	}

	// Evaluate spline i = min at z and return.
	// If we're beyond the points, continue along the tangent
	i = min;
	return y(i) + ((y(i+1) - y(i)) / (x(i+1) - x(i))) * (z-x(i));

}


