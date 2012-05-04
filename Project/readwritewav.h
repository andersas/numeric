#ifndef MYWAVLIB
#define MYWAVLIB
#include <stdint.h>

#define WAVREAD 1
#define WAVWRITE 2

// See the canonical wave file format
// https://ccrma.stanford.edu/courses/422/projects/WaveFormat/

typedef struct {

	char ChunkID[4]; // RIFF
	uint32_t ChunkSize;
	char Format[4]; // Wave

} riff_header;

typedef struct {
	char SubchunkID[4]; // "fmt "
	uint32_t SubchunkSize;
	uint16_t AudioFormat; // Should be 1 for PCM
	uint16_t NumChannels; // 1 = mono, 2 = stereo etc
	uint32_t SampleRate;
	uint32_t ByteRate; // = SampleRate * NumChannels * BitsPerSample/8
	uint16_t BlockAlign; // = NumChannels * BitsPerSample/8
	uint16_t BitsPerSample; // 8, 16, 24 ? ? ?
} riff_fmt_subchunk;

typedef struct {
	char SubchunkID[4]; // "data"
	uint32_t SubchunkSize;

	uint32_t start; // Index of the data, for internal use.
} riff_data_subchunk;


typedef struct {

	riff_header header;
	riff_fmt_subchunk fmt;
	riff_data_subchunk data;
	FILE *fp;
	int mode; // wavread/write

} WAVFILE;


// This code will not work on a big-endian system because
// we read integers directly from a file.
WAVFILE *openwav(char *filename);
WAVFILE *mkwav(char *filename, uint16_t NumChannels, uint32_t SampleRate, uint16_t BitsPerSample);
void wavsync(WAVFILE *);
void wavclose(WAVFILE *F);


// Return 0 on success and -1 on error.
// Data should be an array of pointers to a number of
// arrays for containing the samples.

int64_t readwav(WAVFILE *wavfile, double **data, uint32_t startsample, uint32_t n_samples);
int writewav(WAVFILE *wavfile, double **data, uint32_t startsample, uint32_t n_samples);

// Function to allocate arrays for data to channels.
double **allocate_channels(WAVFILE *F, uint32_t num_samples);
void free_channels(double **pp);


#endif
