#ifndef STEP
#define STEP


double rkstep1(void (*f)(double x, const double y[], double result[]), int dim, double x0, double y0[dim], double h);

double rkstep4(void (*f)(double x, const double y[], double result[]), int dim, double x0, double y0[dim], double h);


#endif
