#ifndef _NEWTON
#define _NEWTON
#include "../linalg/linalg.h"

// Solves f(x) = 0 with newtons method, where f is a function
// from R^N -> R^N
// f is an array with the N coordinate functions of f,
// with x an initial guess. x should be a column vector of length N,
// and dx should be a size that's small compared to the scale of the problem.
// acc is the desired accuracy.

// Returns the number of iterations used, with the root returned in x.
// If this count is greater than 5000, usually something bad has happened
// and the method has failed to converge.

int newton(size_t n, double (*f[n])(matrix *x), matrix *x, double dx, double acc);


#endif
