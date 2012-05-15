#include "simulated_annealing.h"
#include <stdio.h>

double f(double *x) {

	return (*x)*(*x);

}

int main() {

	double x,y;
	double low = 0.5;
	double high = 100.0;

	y = simulated_annealing(1,&x,&f,&low,&high,5.0,1,0.0003,100);

	printf("Min of x^2 = %f, x = %f\n", y, x);

return 0;
}
