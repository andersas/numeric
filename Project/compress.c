#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "readwritewav.h"
#include "readwriteff.h"
#include "compress.h"
#include "../9.FFT/fft.h"

int main(int argc, char **argv) {

double rmax, imax;
double FrequencyResolution;
uint32_t NumChannels, NumSamples, exponent;
uint32_t maxSample, minSample;
uint32_t SampleRate;
uint64_t ReadSamples;
uint32_t i, n, k;
uint32_t max_samples;
double **channels;
complex **channelspreFFT, **channelspostFFT;
FF_data_header header;
WAVFILE *wav;
FFFILE *ff;

if (argc == 1 || argc > 3) {
	printf("Usage: %s in.wav [out.ff]\n", argv[0]);
	printf("Compresses in.wav\n");
	exit(1);
}
if ((wav = openwav(argv[1])) == NULL) {
	printf("Could not open %s...\n", argv[1]);
	exit(2);
}
if (argc == 3) {
	if ((ff = mkff(argv[2], wav->fmt.NumChannels, wav->fmt.SampleRate, wav->fmt.BitsPerSample)) == NULL) {
		printf("Could not open %s...\n", argv[2]);
		exit(3);
	}
} else {
	char filename[strlen(argv[1]) + 3 + 1];
	size_t end;
	strcpy(filename,argv[1]);
	end = strlen(filename);
	// Append .ff to the filename
	filename[end] = '.';
	filename[end+1] = 'f';
	filename[end+2] = 'f';
	filename[end+3] = '\0';
	if ((ff = mkff(filename, wav->fmt.NumChannels, wav->fmt.SampleRate, wav->fmt.BitsPerSample)) == NULL) {
		printf("Could not open %s...\n", filename);
		exit(4);
	}
}

NumChannels = wav->fmt.NumChannels;
SampleRate = wav->fmt.SampleRate;

// Load a number of samples and do a FFT. Throw away
// frequencies greater than 17 kHz and lower than 60 Hz.

// We are dealing with two windows here. The first
// is the time frame wherein we do our FFT, and
// the second is the corrosponding amount of samples needed to do so.
// They are directly related by time = NumSamples/SampleRate
// The window size should be chosen so that
// there is good enough temporal resolution
// to distinguish 1700 kHz from 1600 kHz and 60 Hz from 70 Hz, say.

// The k'th fourier coefficient represents frequency f_ḱ = k/2N,
// where N is the number of samples if k <= N/2
// in the real case and for k > N/2 c[k] = c[N-k]*
// Therefore, 10 Hz should be the inverse of 0.1 seconds = NumSamples/2SampleRate
// NumSamples = 0.1 seconds * 2*SampleRate.
// At 44,100 samples/second (used on most audio CD's) this all translates
// to NumSamples = 4410. To get the nearest power of two,
// we take exponent = (int) log2(NumSamples) + 1 and calculate 2^exponent.
// In this example 8192 samples, with a frequency resolution of 10.77 Hz

exponent = (int) (1.0 + log2(SampleRate/FREQUENCY_RESOLUTION));
if (exponent < 0) {
fprintf(stderr, "Samplerate < 1/2 Frequency resolution.\n");
exit(-1);
}
NumSamples = 1<<exponent; 

channels = allocate_channels(wav, NumSamples);
channelspreFFT = (complex **) malloc(NumChannels*sizeof(complex *));
channelspostFFT = (complex **) malloc(NumChannels*sizeof(complex *));
for (i = 0; i < NumChannels; i++) {
	channelspreFFT[i] = (complex *) malloc(NumSamples*sizeof(complex));
	channelspostFFT[i] = (complex *) malloc(NumSamples*sizeof(complex));
}


FrequencyResolution = 2*SampleRate/NumSamples;
maxSample = MAX_FREQ/FrequencyResolution;
minSample = MIN_FREQ/FrequencyResolution;

// We throw away the upper half of the coefficients (because
// they are mirroring the lower half because we deal with real numbers)
// and the frequencies above MAX_FREQ.
// After throwing away the upper frequencies, in this case we are
// left with 1574 frequency samples.
// Each sample consists of a real part and an imaginary part,
// each BitsPerSample/8 bytes long. Therefore, we have a compression
// rate of 62 % already.
//
// If we then go to throwing away samples of neglicible amplitude,
// we might get even better compression. However, we need a map
// of which samples have been thrown away.
// In this example, we have 1574 samples. If each are represented
// by 1 bit in a bitmap, (1 => kept, 0 => thrown away)
// we need 197 bytes. Each sample takes up 2*BitsPerSample/8 bytes,
// so (worst case: 8 bit/sample) we need to throw away at least
// 12, or 7.8 % of all the remaining samples, before the map takes up
// less space than gained by throwing away samples in the first place.
// This seems resonable in an average sound file.

bitmap_handler_alloc(&header, maxSample-minSample, NumChannels);

max_samples = wav->data.SubchunkSize / (wav->fmt.NumChannels*wav->fmt.BitsPerSample/8);
//printf("%i, %i, %i\n", NumSamples,  minSample, maxSample);
i = 0;
printf("     ");
while ((ReadSamples = readwav(wav,channels,i, NumSamples)) != 0) {

if (i < 0) {
	fprintf(stderr, "Error reading wav file. Aborting\n");
	exit(-1);
}


rmax = imax = 0;
for (n = 0; n < NumChannels; n++) {
	if (ReadSamples < NumSamples) { // Fill the rest with 0's
	for (k = ReadSamples; k < NumSamples; k++)
		channels[n][k] = 0;
	}
	// Prepare for FFT
	for (k = 0; k < NumSamples; k++) {
		channelspreFFT[n][k] = (complex) channels[n][k];
	}

	// FFT:
	fft(NumSamples, channelspreFFT[n], channelspostFFT[n], -1);

	// Readout:
	for (k = minSample; k < maxSample; k++) {
	if (cabs(channelspostFFT[n][k]) > MIN_AMPLITUDE) { // Check if we can throw away sample
		bitmap_handler_set(&header,k-minSample, n);

		channels[n][2*(k-minSample)] = creal(channelspostFFT[n][k]);
		channels[n][2*(k-minSample) + 1] = cimag(channelspostFFT[n][k]);
		
		if (rmax < fabs(channels[n][2*(k-minSample)]))
			rmax = fabs(channels[n][2*(k-minSample)]);
		if (imax < fabs(channels[n][2*(k-minSample) + 1]))
			imax = fabs(channels[n][2*(k-minSample) + 1]);
	} else {
		bitmap_handler_clr(&header,k-minSample, n);
	}
	}
}
header.real_max = rmax;
header.imag_max = imax;
if (ffwritenext(ff, &header, channels, 2*(maxSample-minSample)) != 0) {
fprintf(stderr,"Could not write!\n");
exit(-1);
}
i += ReadSamples;
printf("\b\b\b\b\b");
printf("%3.0f %%", ((double) i) / ((double) max_samples) * 100.0);
fflush(stdout);
}
printf("\n");
wavclose(wav);
wavclose(ff);


return 0;
}



