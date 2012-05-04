#ifndef FFLIB
#define FFLIB
#include <stdint.h>
#include "readwritewav.h"

#define FFREAD 1
#define FFWRITE 2

typedef WAVFILE FFFILE;

typedef struct {
	float real_max;
	float imag_max;
	int num_channels;
	uint32_t bitmap_size;
	uint8_t **bitmap;
} FF_data_header;

void bitmap_handler_alloc(FF_data_header *h, uint32_t num_bits, int num_channels);
void bitmap_handler_free(FF_data_header *h);
void bitmap_handler_set(FF_data_header *h, uint32_t bitno, int channel);
void bitmap_handler_clr(FF_data_header *h, uint32_t bitno, int channel);
int bitmap_handler_get(FF_data_header *h, uint32_t bitno, int channel);

// This code will not work on a big-endian system because
// we read integers directly from a file.
void ffrewind(FFFILE *F);
FFFILE *openff(char *filename);
FFFILE *mkff(char *filename, uint16_t NumChannels, uint32_t SampleRate, uint16_t BitsPerSample);
void ffsync(FFFILE *);
void ffclose(FFFILE *F);


// Return 0 on success and -1 on error.
// Data should be an array of pointers to a number of
// arrays for containing the samples.

int ffreadnext(FFFILE *fffile, FF_data_header *header, double **data, uint32_t n_samples);
int ffwritenext(FFFILE *fffile, FF_data_header *header, double **data, uint32_t n_samples);

#endif
