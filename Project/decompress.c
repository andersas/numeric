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

double FrequencyResolution;
uint32_t NumChannels, NumSamples, exponent;
uint32_t maxSample, minSample;
uint32_t SampleRate;
uint64_t ReadSamples;
uint32_t i, n, k;
double **channels;
complex **channelspreFFT, **channelspostFFT;
FF_data_header header;
WAVFILE *wav;
FFFILE *ff;

if (argc == 1 || argc > 3) {
	printf("Usage: %s in.ff [out.wav]\n", argv[0]);
	printf("Decompresses in.ff\n");
	exit(1);
}
if ((ff = openff(argv[1])) == NULL) {
	printf("Could not open %s...\n", argv[1]);
	exit(2);
}
if (argc == 3) {
	if ((wav = mkwav(argv[2], ff->fmt.NumChannels, ff->fmt.SampleRate, ff->fmt.BitsPerSample)) == NULL) {
		printf("Could not open %s...\n", argv[2]);
		exit(3);
	}
} else {
	char filename[strlen(argv[1]) + 4 + 1];
	size_t end;
	strcpy(filename,argv[1]);
	end = strlen(filename);
	// Append .wav to the filename
	filename[end] = '.';
	filename[end+1] = 'w';
	filename[end+2] = 'a';
	filename[end+3] = 'v';
	filename[end+4] = '\0';
	if ((wav = mkwav(filename, ff->fmt.NumChannels, ff->fmt.SampleRate, ff->fmt.BitsPerSample)) == NULL) {
		printf("Could not open %s...\n", filename);
		exit(4);
	}
}

NumChannels = ff->fmt.NumChannels;
SampleRate = ff->fmt.SampleRate;

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

exponent = 1.0 + log2(SampleRate/FREQUENCY_RESOLUTION);
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
// they are mirroring the lower half) and the frequencies
// above MAX_FREQ.
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
// 197, or 7.8 % of all the remaining samples, before the map takes up
// less space than gained by throwing away samples in the first place.
// This seems resonable in an average sound file.

bitmap_handler_alloc(&header, maxSample-minSample, NumChannels);

i = 0;
while ((ReadSamples = ffreadnext(ff,&header,channels,2*(maxSample-minSample))) != 0) {
if (ReadSamples < 0) {
	fprintf(stderr, "Error reading wav file. Aborting\n");
	exit(-1);
}

for (n = 0; n < NumChannels; n++) {
	 for (k = 0; k <= NumSamples/2; k++) {
                if (k < minSample || k >= maxSample) {
                        channelspreFFT[n][k] = 0;
                } else {
                channelspreFFT[n][k] = channels[n][2*(k-minSample)] + I * channels[n][2*(k-minSample) + 1];
                }
        }
        for (k = NumSamples/2+1; k < NumSamples; k++) {
                channelspreFFT[n][k] = conj(channelspreFFT[n][NumSamples-k]);
        }

        fft(NumSamples, channelspreFFT[n], channelspostFFT[n], 1);

        for (k = 0; k < NumSamples; k++) {
                channels[n][k] = creal(channelspostFFT[n][k]);
        }
}
if (writewav(wav, channels, i, NumSamples) != 0){
fprintf(stderr,"Could not write!\n");
exit(-1);
}

i += NumSamples;
}
printf("\n");
wavclose(wav);
wavclose(ff);


return 0;
}



