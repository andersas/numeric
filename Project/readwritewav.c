#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "readwritewav.h"

#define min(a,b) (a < b) ? a : b;

// The file must be about to be read just after
// the FMT header.

int parse_fmt(WAVFILE *F) { 


	if (fread(&(F->fmt.AudioFormat), 2, 1, F->fp) != 1) return 1;
	if (fread(&(F->fmt.NumChannels), 2, 1, F->fp) != 1) return 1;
	if (fread(&(F->fmt.SampleRate), 4, 1, F->fp) != 1) return 1;
	if (fread(&(F->fmt.ByteRate), 4, 1, F->fp) != 1) return 1;
	if (fread(&(F->fmt.BlockAlign), 2, 1, F->fp) != 1) return 1;
	if (fread(&(F->fmt.BitsPerSample), 2, 1, F->fp) != 1) return 1;

	if (F->fmt.AudioFormat != 1) return 1; // PCM
	if (F->fmt.NumChannels == 0) return 1; // No channels
	if (F->fmt.ByteRate != F->fmt.SampleRate * F->fmt.NumChannels * F->fmt.BitsPerSample / 8) return 1;
	if (F->fmt.BlockAlign != F->fmt.NumChannels * F->fmt.BitsPerSample / 8) return 1;

	if (8 * (F->fmt.BitsPerSample / 8) != F->fmt.BitsPerSample || F->fmt.BitsPerSample > 32) return 1;

	return 0;
}

WAVFILE *openwav(char *filename) {

	// This code will not work on a big-endian system because
	// we read integers directly from a file.
	
	struct stat statbuf;
	WAVFILE *wavfile;
	FILE *file;
	char ChunkID[4];
	size_t size;
	size_t position;
	uint32_t next_size;
	int num_fmt = 0; // Counts number of subchunks
	int num_data = 0;

	if ((file = fopen(filename, "r")) == NULL)
		return NULL;

	wavfile = (WAVFILE *) malloc(sizeof(WAVFILE));

	wavfile->mode = WAVREAD;	
	wavfile->fp = file;	

	stat(filename, &statbuf);
	size = statbuf.st_size;
	if (size < 44) goto error; // Minimum size for empty wav file with header
	// Parse the RIFF header. Should start with RIFF<size>WAVE
	if (fread(&(wavfile->header.ChunkID), 4, 1, file) != 1) goto error;
	if (fread(&(wavfile->header.ChunkSize), 4, 1, file) != 1) goto error;
	if (fread(&(wavfile->header.Format), 4,1,file) != 1) goto error;
	if (strncmp(wavfile->header.ChunkID, "RIFF", 4)) goto error;
	if (strncmp(wavfile->header.Format, "WAVE", 4)) goto error;

	if (wavfile->header.ChunkSize % 2 == 0) {
		if (wavfile->header.ChunkSize != size - 8) {
			 goto error;
		}
	} else if (wavfile->header.ChunkSize + 1 != size - 8)  {
		goto error;
	}
	
	// Parse the subchunks of the wavefile
	position = 12;
	
	do {
	if (fread(ChunkID, 4, 1, file) != 1) goto error;
	if (fread(&next_size, 4, 1, file) != 1) goto error;
//printf("%i: \"%c%c%c%c\": %i\n", position, ChunkID[0], ChunkID[1], ChunkID[2], ChunkID[3], next_size);
		if (next_size + position + 4 > size) goto error;
		if (!strncmp(ChunkID,"fmt ", 4)) {
			if (num_fmt == 0 && next_size >= 16) {
				strncpy(wavfile->fmt.SubchunkID, ChunkID, 4);
				wavfile->fmt.SubchunkSize = next_size;
				if (parse_fmt(wavfile)) goto error;
				num_fmt++;
			} else
				goto error;
		} else if (!strncmp(ChunkID, "data", 4)) {
			if (num_data == 0) {
                                strncpy(wavfile->data.SubchunkID, ChunkID, 4);
                                wavfile->data.SubchunkSize = next_size;
				wavfile->data.start = position + 8;
                                num_data++;
                        } else
				goto error;
		}
	position += next_size + 8; // each header is 8 bytes long
	if (position % 2 == 1) position++;
	fseek(file, (long) position, SEEK_SET);
	} while (position < size && position > 16); // the last check
						    // is to see if we've
						    // wrapped around.

	return wavfile;

error:
	free(wavfile);
	fclose(file);
	return NULL;
}



WAVFILE *mkwav(char *filename, uint16_t NumChannels, uint32_t SampleRate, uint16_t BitsPerSample) {

	WAVFILE *wavfile;
	FILE *file;

	if (8 * (BitsPerSample / 8) != BitsPerSample || BitsPerSample > 32) return NULL;

	if ((file = fopen(filename, "w")) == NULL) {
		fprintf(stderr,"%s\n", filename);
		perror("fopen");
		return NULL; 
	}

	wavfile = (WAVFILE *) malloc(sizeof(WAVFILE));

	wavfile->mode = WAVWRITE;
	wavfile->fp = file;

	strncpy(wavfile->header.ChunkID, "RIFF", 4);
	strncpy(wavfile->header.Format, "WAVE", 4);
	strncpy(wavfile->fmt.SubchunkID, "fmt ", 4);
	strncpy(wavfile->data.SubchunkID, "data", 4);
	
	wavfile->header.ChunkSize = 24+8; // No data in data chunk
	wavfile->fmt.SubchunkSize = 24-8;
	wavfile->data.SubchunkSize = 0;
	
	wavfile->fmt.AudioFormat = 1; // PCM
	wavfile->fmt.NumChannels = NumChannels;
	wavfile->fmt.SampleRate = SampleRate;
	wavfile->fmt.ByteRate = SampleRate * NumChannels * BitsPerSample/8;
	wavfile->fmt.BlockAlign = NumChannels * BitsPerSample/8;
	wavfile->fmt.BitsPerSample = BitsPerSample;
	wavfile->data.start = 44;

	wavsync(wavfile);

	return wavfile;
}

// Writes the header to disk
void wavsync(WAVFILE *wavfile) {


	if (wavfile->mode != WAVWRITE) {
		fprintf(stderr, "Trying to write to readonly wavfile.\n");
		return;
	}
	fseek(wavfile->fp, 0, SEEK_SET);
	if (fwrite(wavfile->header.ChunkID, 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->header.ChunkSize), 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(wavfile->header.Format, 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(wavfile->fmt.SubchunkID, 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->fmt.SubchunkSize), 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->fmt.AudioFormat), 2, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->fmt.NumChannels), 2, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->fmt.SampleRate), 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->fmt.ByteRate), 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->fmt.BlockAlign), 2, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->fmt.BitsPerSample), 2, 1, wavfile->fp) != 1) goto error;
	if (fwrite(wavfile->data.SubchunkID, 4, 1, wavfile->fp) != 1) goto error;
	if (fwrite(&(wavfile->data.SubchunkSize), 4, 1, wavfile->fp) != 1) goto error;

	return;

error:
	fprintf(stderr, "Wavsync could not write to open wavfile.\n");
	exit(-1);

}

void wavclose(WAVFILE *F) {

	fclose(F->fp);
	free(F);

}

// Return number of success samples read and -1 on error.
// Data should be an array of pointers to a number of
// arrays for containing the samples.

int64_t readwav(WAVFILE *wavfile, double **data, uint32_t startsample, uint32_t n_samples) {

	FILE *file = wavfile->fp;
	uint16_t NumChannels = wavfile->fmt.NumChannels;
	uint16_t BitsPerSample = wavfile->fmt.BitsPerSample;
	uint32_t start = wavfile->data.start;
	void *block = malloc(NumChannels*(BitsPerSample/8));
	uint32_t i,n;
	int32_t int24t;
	uint32_t readsamples = 0;

	n_samples = min(n_samples, wavfile->data.SubchunkSize / (NumChannels*(BitsPerSample/8)));

	fseek(file, (long) start + startsample*(BitsPerSample/8)*NumChannels, SEEK_SET);

	for (i = 0; i < n_samples; i++) {

	if (wavfile->data.start + (i+startsample)*(BitsPerSample/8)*NumChannels > wavfile->data.SubchunkSize) break;

		if (fread(block, (BitsPerSample/8)*NumChannels, 1, file) != 1) {
			break;
		}

		for (n = 0; n < NumChannels; n++) {
			if (BitsPerSample == 8) {
				data[n][i] = ((double) ((int8_t *) block)[n] - 128) / 127;
			} else if (BitsPerSample == 16) {
				data[n][i] = ((double) ((int16_t *) block)[n]) / 32768.0; // 2^15
			} else if (BitsPerSample == 24) {
				int24t = (int32_t) ( \
				    (((uint8_t *) block)[3*n]) + \
				    ((((uint8_t *) block)[3*n + 1]) << 8) + \
				    ((((uint8_t *) block)[3*n + 2]) << 16));
				// Sign extend
				data[n][i] = (((double) int24t) * \
					((int24t & 0x00800000) ? -1 : 1));
				data[n][i] /= 8388608.0; //2^23 
					       
			} else if (BitsPerSample == 32) {
				data[n][i] = ((double) ((int32_t *) block)[n]) / 2147483648.0; // 2^31
			} else {
				free(block);
				return -1;
			}

		}
		readsamples = i + 1;
	}

	free(block);
	return readsamples;
}

// Return 0 on success and -1 on error.
// Data should be an array of pointers to a number of
// arrays for containing the samples.

// Write to a wavfile assumed to be opened by mkwav.

int writewav(WAVFILE *wavfile, double **data, uint32_t startsample, uint32_t n_samples) {

// Again, we assume we're on a little-endian system (like the x86 family of
// processors)

	FILE *file = wavfile->fp;
	uint16_t NumChannels = wavfile->fmt.NumChannels;
	uint16_t BitsPerSample = wavfile->fmt.BitsPerSample;
	uint32_t start = wavfile->data.start;
	void *block = malloc(NumChannels*(BitsPerSample/8));
	uint32_t i,n;
	int32_t int24t;
	int64_t diff;

	if (wavfile->mode != WAVWRITE) {
		fprintf(stderr, "Attempting to write read-only wav file\n");
		return -1;
	}

	diff = startsample + n_samples - wavfile->data.SubchunkSize / (NumChannels*(BitsPerSample/8));

	if (diff > 0) {
		// Update headers to reflect new size 
		wavfile->data.SubchunkSize = (NumChannels*(BitsPerSample/8)) * (startsample + n_samples);
		wavfile->header.ChunkSize = 4 + 24 + 8 + wavfile->data.SubchunkSize;
		wavfile->fmt.BlockAlign = NumChannels*(BitsPerSample/8);
		if (ftruncate(fileno(wavfile->fp), 8 + wavfile->header.ChunkSize)) {
			perror("ftruncate");
			exit(-1);
		}
		wavsync(wavfile); // write actual headers to disk
	}


	fseek(file, (long) start + startsample*(BitsPerSample/8)*NumChannels, SEEK_SET);

	for (i = 0; i < n_samples; i++) {
		for (n = 0; n < NumChannels; n++) {
		if (BitsPerSample == 8) {
			((uint8_t *) block)[n] = (uint8_t) (data[n][i]*127 + 128); 
		} else if (BitsPerSample == 16) {
			((uint16_t *) block)[n] = (uint16_t) (data[n][i]*32767.0);
		} else if (BitsPerSample == 24) {
			int24t = (int32_t) data[n][i] * 8388607.0;
			((uint8_t *) block)[n] = (uint8_t) (int24t & 0xFF);
			((uint8_t *) block)[n+1] = (uint8_t) (int24t & 0xFF00) >> 8;
			((uint8_t *) block)[n+2] = (uint8_t) (int24t & 0xFF0000) >> 16;
		} else if (BitsPerSample == 32) {
			((uint32_t *) block)[n] = (uint32_t) (data[n][i]*2147483647.0);
		}

		}

		if (fwrite(block, NumChannels*(BitsPerSample/8), 1, file) != 1) {
			fprintf(stderr,"Could not write to wavfile. It might be corrupt.\n");
			exit(-1);
		}
	}

	free(block);
	return 0;
}




// Function to allocate arrays for data to channels.
double ** allocate_channels(WAVFILE *F, uint32_t num_samples) {

	uint32_t i;
	double **pp;
	double **pp2;

	pp = (double **) malloc(sizeof(double) + F->fmt.NumChannels * sizeof(double *));
	pp2 = &pp[1];

	for (i = 0; i < F->fmt.NumChannels; i++) {
		pp2[i] = (double *) malloc(num_samples*sizeof(double));
	}
	
	// This is somewhat of a hack to get a header in
	pp[0] = (double *) malloc(sizeof(uint32_t)); 
	memcpy(pp[0], &num_samples, sizeof(num_samples));

	return pp2;
}

void free_channels(double **pp) {
	uint32_t num_samples;
	uint32_t i;

	double **pp2 = &pp[1];
	memcpy(&num_samples, pp[0], sizeof(num_samples));

	for (i = 0; i < num_samples; i++) {
		free(pp2[i]);
	}
	free(pp[0]);

	free(pp);
}
