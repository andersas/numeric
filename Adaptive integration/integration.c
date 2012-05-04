#include <stdio.h>
#include <math.h>
#include "adaptive.h"

double f1(double x) {
	return log(x)/sqrt(x);
}
double f2(double x) {
	return 4*sqrt(1.-(1.-x)*(1.-x));
}

double step(double x) {
	return x > 0 ? 1.0 : 0.0;
}

double fness(double x) {

	return ( 1+0.2*sin ( x ) ) * ( x*x+1 ) *exp ( -x*x/16 );

}

int main () {

	struct NC24result integration;

	printf("Integrating log(x)/sqrt(x) from 0 to 1...\n");

	integration = NC24integrate(&f1,0,1,0.001, 0.001);

	printf("Result: %lf, error: %lf\n", integration.result, integration.error);

	printf("\n");

	printf("Integrating 4sqrt(1-(1-x)^2) from 0 to 1...\n");

	integration = NC24integrate(&f2,0,1,1e-11,0.0); 


	printf("M_PI =  %.25f\n",M_PI);
	printf("pi   =  3.1415926535897932384626433832795028841971693993751058209749445923078164\n");
	printf("Result: %.17lf,\nerror:  %.17lf\n", integration.result, integration.error);


	printf("\nIntegrating step function from -1 to 1...\n");
	integration = NC24integrate(&step,-1,1,0.001, 0.001);
	printf("Result: %lf, error: %lf\n", integration.result, integration.error);

	printf("\nIntegrating weird function from -10 to 10...\n");
	integration = NC24integrate(&fness,-10,10,0.001, 0.001);
	printf("Result: %lf, error: %lf\n", integration.result, integration.error);


	return 0;


}
