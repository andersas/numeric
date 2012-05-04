#include <stdio.h>
#include <complex.h>
#include <string.h>
#include "fft.h"

int main(int argc, char **argv) {

	complex x[4096];
	complex c[4096];
	double buf;
	char *filename;
	int do_dft = 0;
	int i;
	FILE *inp, *outp;


	if (argc > 1) {
	if (!strcmp(argv[1], "dft")) {
	do_dft = 1;
	}}
	if (do_dft)
	printf("Reading 4096 points from inp.dat and perfoming dft\n");
	else
	printf("Reading 4096 points from inp.dat and perfoming fft\n");

	filename = "out-fft.dat";
	if (do_dft)
		filename = "out-dft.dat";

	inp = fopen("inp.dat", "r");
	outp = fopen(filename, "w");

	for ( i = 0; i < 4096; i++) {
		fscanf(inp,"%lf\n", &buf);
		x[i] = (complex) buf;
	}

	if (!do_dft)
		fft(4096,x,c,-1);
	else 
		dft(4096,x,c,-1);


	printf("Dumping 4096 transformed points to %s...\n",filename);
	for ( i = 0; i < 4096; i++) {
		fprintf(outp,"%i %f %fi\n", i, creal(c[i]), cimag(c[i]));
	}

	fclose(inp);
	fclose(outp);

	if (argc == 3) {
		filename = "out-inverse-fft.dat";
		if (do_dft)
			filename = "out-inverse-dft";
		printf("Doing invers transform, dumping to %s.\n",filename);

		if (do_dft)
			dft(4096, c,x,1);
		else
			fft(4096,c,x,1);

	outp = fopen(filename, "w");
	for (i = 0; i < 4096; i++)
		fprintf(outp,"%i %f %fi\n", i, creal(x[i]), cimag(x[i]));

	fclose(outp);
	}

return 0;
}

