#include <stdio.h>
#include <stdlib.h>
#include "lspline.h"
#include "data.h"


#define ALLOC_INCREMENT 256


int main(int argc, char **argv) {

double x,y,start,stop,step;
int n = 0;

struct point * input;

if (argc != 4) {
	fprintf(stderr, "Usage: %s start stop step\n\n", argv[0]);
	fprintf(stderr, "Prints interpolated function points from start to stop\n");
	fprintf(stderr, "in steps of size step of input data on stdin.\n");
	exit(2);
}

start = atof(argv[1]);
stop = atof(argv[2]);
step = atof(argv[3]);

if (stop < start) {
	fprintf(stderr, "stop < start\n");
	exit(3);
}
if (step < 0.00000001) {
	fprintf(stderr, "Too small step size\n");
	exit(4);
}

input = (struct point *) malloc(sizeof(struct point) * ALLOC_INCREMENT);

/* Read input data from standard input */
while (scanf("%lf %lf", &x, &y) != EOF) {

	if (n / ALLOC_INCREMENT > (n-1)/ALLOC_INCREMENT) {
		/* Allocate more space, as we've run out */
		input = (struct point *) realloc((void *) input, \
	 		 sizeof(struct point) * (n + ALLOC_INCREMENT));
	}

	input[n].x = x;
	input[n].y = y;
	n++;
}


/* Then perform linear spline interpolation and print the result */

while (start <= stop) {
	printf("%lf\t%lf\n", start, lspline_interpolate(n,input,start));
	start += step;
}

return 0;


}
