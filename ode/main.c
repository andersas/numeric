#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "driver.h"

int count[5];
int count1[5];
int count4[5];

void example1(double x, const double y[], double result[]) {
	// Basic example: sine (should give us -cos(x))
	count[0]++;
	result[0] = sin(x);
}
#define example1_start 0.0
#define example1_stop 10.0
double example1_y0[] = {-1.0};
#define example1_dim 1

void example2(double x, const double y[], double result[]) {
	// Coupled differential equations
	// y1' = 10(y2-y1)
	// y2' = y1(28-y3)-y2
	// y3' = y1y2 - 8/3y3
	count[1]++;
	result[0] = 10*y[1]-10*y[0];
	result[1] = y[0]*28 - y[0]*y[2] - y[1];
	result[2] = y[0]*y[1] - 8.0*y[2]/3.0;
	
}
#define example2_start 0.0
#define example2_stop 100.0
double example2_y0[] = {3.1,7.9,1.8};
#define example2_dim 3

void example3(double x, const double y[], double result[]) {
	// Logistic growth
	count[2]++;
	result[0] = y[0]*(1-y[0]);
}
#define example3_start 0.0
#define example3_stop 20.0
double example3_y0[] = {0.00001};
#define example3_dim 1

void example4(double x, const double y[], double result[]) {
	// Second order ODE (dampened harmonic motion)
	// x'' = -cx' - x
	// Setting
	// y0 = x,
	// y1 = x'
	// y2 = x''
	// we get
	// y0' = y1
	// y1' = y2
	// y2' = -c*y2 - 1
	
	double criticallity = 0.1;// > 1 over dampening, = 1 critical, < 1 under
	count[3]++;
	result[0] = y[1];
	result[1] = y[2];
	result[2] = -criticallity * y[2] - y[1];
}
#define example4_start 0.0
#define example4_stop 100.0
double example4_y0[] = {10.0, 10.0, 10.0};
#define example4_dim 3


void example5(double x, const double y[], double result[]) {
	// Something like a random walk
	count[4]++;
	result[0] = (random() > 1073741823) - 0.5;
}
#define example5_start 0.0
#define example5_stop 10.0
double example5_y0[] = {0.0};
#define example5_dim 1

double start[] = {
	example1_start,
	example2_start,
	example3_start,
	example4_start,
	example5_start
};
double stop[] = {
	example1_stop,
	example2_stop,
	example3_stop,
	example4_stop,
	example5_stop
};
double *init[] = {
	example1_y0,
	example2_y0,
	example3_y0,
	example4_y0,
	example5_y0
};
int dim[] = {
	example1_dim,
	example2_dim,
	example3_dim,
	example4_dim,
	example5_dim
};

void (*example[5])(double x, const double y[], double result[]) = {
	&example1,
	&example2,
	&example3,
	&example4,
	&example5
};

FILE *file1[5]; // Return values
FILE *file4[5]; // Return values


// Run an example with rk1 and rk4 and print it out to a file,
// reporting how many function calls it took each time.
void * run_example(void *num) {

	const int i = (int) num;
	int n;
	double buf[dim[i]];
	driver_r *r, *tmp;

	count[i] = 0;
	memcpy(buf, init[i], dim[i]*sizeof(double));
	tmp = r = driver(example[i], dim[i], start[i], stop[i], init[i], 0.01, 0.01, 0.1, 1);
	memcpy(init[i], buf, dim[i]*sizeof(double));
	while (r != NULL) {
		fprintf(file1[i], "%f", r->x);
		for (n = 0; n < dim[i]; n++)
			fprintf(file1[i], "\t%f", r->y[n]);
		fprintf(file1[i], "\n");
		r = r->next;
	}
	count1[i] = count[i];
	driver_r_destroy(tmp);


	count[i] = 0;
	tmp = r = driver(example[i], dim[i], start[i], stop[i], init[i], 0.01, 0.01, 0.1, 4);
	while (r != NULL) {
		fprintf(file4[i], "%f", r->x);
		for (n = 0; n < dim[i]; n++)
			fprintf(file4[i], "\t%f", r->y[n]);
		fprintf(file4[i], "\n");
		r = r->next;
	}
	count4[i] = count[i];
	driver_r_destroy(tmp);

	return NULL;

}

int main(int argc, char **argv) {

	char path[20];
	pthread_t thread[5];
	int i;
	pthread_attr_t attr;

// Open files for writing
	strcpy(path, "example");
	for (i = 0; i < 5; i++) {
		sprintf(&path[7], "-%i-rk1.dat", i+1);
		file1[i] = fopen(path, "w");
		sprintf(&path[7], "-%i-rk4.dat", i+1);
		file4[i] = fopen(path, "w");
	}
// Run examples in parallel

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for (i = 0; i < 5; i++) {
		pthread_create(&thread[i], &attr, &run_example, (void *) i);
	}

	pthread_attr_destroy(&attr);

// Close and exit

	for (i = 0; i < 5; i++) {
		pthread_join(thread[i], NULL); 
		fclose(file1[i]);
		fclose(file4[i]);
		printf("Example %i done with rk1 took %i calls.\n", i + 1, count1[i]);
		printf("Example %i done with rk4 took %i calls.\n", i + 1, count4[i]);
	}


return 0;

}
