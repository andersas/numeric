#include <stdio.h>
#include <math.h>
#include "adaptive.h"

// Estimate the integral of f(x) from a to b to within
// an absolute/relative error of acc/eps
// using the recursive adaptive open 2/4 Newton-Cotes quadrature 

int calls, depth, maxdepth;

struct NC24result integrate(double (*f)(double x), double a, double b, double acc, double eps, double save1, double save2) {

	struct NC24result result;
	double fs[4]; // Saved f values
	double h = b - a;
	double q4, q2;

//	Abscissas for 4th order quadrature:
//	(refer to eq. 9 and 10 in lecture notes)
//	{1.0/6.0, 2.0/6.0, 4.0/6.0, 5.0/6.0};
//	Bascissas for 2nd order quadrature:
//	{1/3, 2/3} (= {2/6, 4/6} happily)

	calls++;  // Count recursion depth etc.
	depth++;
	if (depth > maxdepth)
		maxdepth = depth;

	// We are given f at the inner points
	// (a + h*1/3 and a + h*2/3).
	// Calculate the outer points 

	fs[0] = f(a + 1./6.*h);
	fs[1] = save1;
	fs[2] = save2;
	fs[3] = f(a + 5./6.*h);

	// Evaluate quadratures
	
	q4 = h*(2*fs[0] + fs[1] + fs[2] + 2*fs[3])/6.0;
	q2 = h*(fs[1] + fs[2])/2.0;

	// See if we get an acceptable error
	// (see the lecture notes)
	if (fabs(q4-q2) < acc + eps*fabs(q4)) {
		result.result = q4;
		result.error = fabs(q4 - q2);
	} else { // We didn't, so evaluate integral recursively.
		struct NC24result left, right;
		// absolute error should be divided by sqrt(2)
		// since we ask for half of the integral.
		left  = integrate(f,a,(a+b)/2,acc/M_SQRT2,eps, fs[0], fs[1]);
		right = integrate(f,(a+b)/2,b,acc/M_SQRT2,eps, fs[2], fs[3]);

		result.result = left.result + right.result;
		result.error = sqrt(left.error*left.error + right.error*right.error);
	}

	depth--;
	return result;
}

struct NC24result NC24integrate(double (*f)(double x), double a, double b, double acc, double eps)
{

	struct NC24result r;
	double h = b - a;
	maxdepth = depth = calls = 0;
	r = integrate(f,a,b,acc,eps,f(a + 2.0/6.0*h), f(a+4.0/6.0*h));

	printf("Maximum recursion depth: %i calls.\n", maxdepth);
	printf("Number of calls: %i.\n", calls);

	return r;
}

