#ifndef LEASTSQUARES
#define LEASTSQUARES
#include "../linalg/linalg.h"

void leastSquares(size_t n, size_t m, matrix *x, matrix *y, matrix *dy, double (*f[m])(double), matrix *coeff, matrix *cov);



#endif
