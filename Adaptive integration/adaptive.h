#ifndef ADAPTIVE
#define ADAPTIVE

struct NC24result {
	double result;
	double error;
};

struct NC24result NC24integrate(double (*f)(double x), double a, double b, double acc, double eps);

#endif
