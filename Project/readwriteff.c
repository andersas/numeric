#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "readwritewav.h"
#include "readwriteff.h"

// Functions for handling .ff files, which store compressed wave files.
// A .ff file is a RIFF file like WAVE with the same FMT block
// but the data block consists of frames, each beginning with a header
// followed by the compressed samples, which are Fourier coefficients.

int parse_fmt(FFFILE *F);

#define min(a,b) (a < b) ? a : b;


FFFILE *openff(char *filename) {

	// This code will not work on a big-endian system because
	// we read integers directly from a file.
	
	struct stat statbuf;
	FFFILE *fffile;
	FILE *file;
	char ChunkID[4];
	size_t size;
	size_t position;
	uint32_t next_size;
	int num_fmt = 0; // Counts number of subchunks
	int num_data = 0;

	if ((file = fopen(filename, "r")) == NULL)
		return NULL;

	fffile = (FFFILE *) malloc(sizeof(FFFILE));

	fffile->mode = FFREAD;	
	fffile->fp = file;	

	stat(filename, &statbuf);
	size = statbuf.st_size;

	if (size < 44) goto error; // Minimum size for empty ff file with header
	// Parse the RIFF header. Should start with RIFF<size>FFE
	if (fread(&(fffile->header.ChunkID), 4, 1, file) != 1) goto error;
	if (fread(&(fffile->header.ChunkSize), 4, 1, file) != 1) goto error;
	if (fread(&(fffile->header.Format), 4,1,file) != 1) goto error;
	if (strncmp(fffile->header.ChunkID, "RIFF", 4)) goto error;
	if (strncmp(fffile->header.Format, "FFTC", 4)) goto error;

	if (fffile->header.ChunkSize % 2 == 0) {
		if (fffile->header.ChunkSize != size - 8) {
			 goto error;
		}
	} else if (fffile->header.ChunkSize + 1 != size - 8)  {
		goto error;
	}
	
	// Parse the subchunks of the ffefile
	position = 12;
	
	do {
	if (fread(ChunkID, 4, 1, file) != 1) goto error;
	if (fread(&next_size, 4, 1, file) != 1) goto error;
//printf("%i: \"%c%c%c%c\": %i\n", position, ChunkID[0], ChunkID[1], ChunkID[2], ChunkID[3], next_size);
		if (next_size + position + 4 > size) goto error;
		if (!strncmp(ChunkID,"fmt ", 4)) {
			if (num_fmt == 0 && next_size >= 16) {
				strncpy(fffile->fmt.SubchunkID, ChunkID, 4);
				fffile->fmt.SubchunkSize = next_size;
				if (parse_fmt(fffile)) goto error;
				num_fmt++;
			} else
				goto error;
		} else if (!strncmp(ChunkID, "data", 4)) {
			if (num_data == 0) {
                                strncpy(fffile->data.SubchunkID, ChunkID, 4);
                                fffile->data.SubchunkSize = next_size;
				fffile->data.start = position + 8;
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
	ffrewind(fffile);
	return fffile;

error:
	free(fffile);
	fclose(file);
	return NULL;
}

void ffrewind(FFFILE *F) {
	fseek(F->fp, F->data.start, SEEK_SET);
}


FFFILE *mkff(char *filename, uint16_t NumChannels, uint32_t SampleRate, uint16_t BitsPerSample) {

	FFFILE *fffile;
	FILE *file;

	if (8 * (BitsPerSample / 8) != BitsPerSample || BitsPerSample > 32) return NULL;

	if ((file = fopen(filename, "w")) == NULL) {
		fprintf(stderr,"%s\n", filename);
		perror("fopen");
		return NULL; 
	}

	fffile = (FFFILE *) malloc(sizeof(FFFILE));

	fffile->mode = FFWRITE;
	fffile->fp = file;

	strncpy(fffile->header.ChunkID, "RIFF", 4);
	strncpy(fffile->header.Format, "FFTC", 4);
	strncpy(fffile->fmt.SubchunkID, "fmt ", 4);
	strncpy(fffile->data.SubchunkID, "data", 4);
	
	fffile->header.ChunkSize = 24+8; // No data in data chunk
	fffile->fmt.SubchunkSize = 24-8;
	fffile->data.SubchunkSize = 0;
	
	fffile->fmt.AudioFormat = 1; // PCM
	fffile->fmt.NumChannels = NumChannels;
	fffile->fmt.SampleRate = SampleRate;
	fffile->fmt.ByteRate = SampleRate * NumChannels * BitsPerSample/8;
	fffile->fmt.BlockAlign = NumChannels * BitsPerSample/8;
	fffile->fmt.BitsPerSample = BitsPerSample;
	fffile->data.start = 44;

	wavsync(fffile);

	return fffile;
}

void ffclose(FFFILE *F) {

	fclose(F->fp);
	free(F);

}


/* Be careful here. n_samples is the number of real valued samples,
 * so each actual sample (a complex number) consists of two samples
 * to be read here... */
int ffreadnext(FFFILE *fffile, FF_data_header *header, double **data, uint32_t n_samples) {

	FILE *file = fffile->fp;
	uint16_t NumChannels = fffile->fmt.NumChannels;
	uint16_t BitsPerSample = fffile->fmt.BitsPerSample;
	uint32_t i,n,m;
	uint32_t skipped_samples;
	int32_t int24t;
	char tmp;
	float max;
	uint32_t readsamples = 0;
	void *block = malloc(NumChannels*(BitsPerSample/8));

// Parse header. 42, frame max real, max imag, bitmap of samples left out.
	if (fread(&tmp, 1,1,file) != 1) {
		return 0;
	}
	if (tmp != 42) {
		fprintf(stderr,"Readnext failed: not aligned.\n");
		exit(-1);
	}

	if (fread(&(header->real_max), sizeof(float),1,file) != 1) return -1;
	if (fread(&(header->imag_max), sizeof(float),1,file) != 1) return -1;
	for (n = 0; n < NumChannels; n++)
		if (fread(header->bitmap[n], header->bitmap_size, 1, file) != 1) return -1;

	for (i = 0; i < n_samples; i++) {
		max = ((i%2 == 0) ? header->real_max : header->imag_max);

		skipped_samples = 0;
		for (n = 0 ; n < NumChannels; n++) {
			if (!bitmap_handler_get(header,i/2,n)) skipped_samples++;
		}
		if (skipped_samples != NumChannels) {
		if (fread(block, (NumChannels-skipped_samples)*(BitsPerSample/8), 1, file) != 1) {
			return -1;
		}
		}
	
		// Data block read in, now parse it.
		skipped_samples = 0;
		for (n = 0; n < NumChannels; n++) {
			if (bitmap_handler_get(header,i/2,n)) {
			m = n - skipped_samples;

			if (BitsPerSample == 8) {
				data[n][i] = ((double) ((int8_t *) block)[m]) / 127.0;
			} else if (BitsPerSample == 16) {
				data[n][i] = ((double) ((int16_t *) block)[m]) / 32767.0; // 2^15 - 1;
			} else if (BitsPerSample == 24) {
				int24t = (int32_t) ( \
				    (((uint8_t *) block)[3*m]) + \
				    ((((uint8_t *) block)[3*m + 1]) << 8) + \
				    ((((uint8_t *) block)[3*m + 2]) << 16));
				// Sign extend
				data[n][i] = (((double) int24t) * \
					((int24t & 0x00800000) ? -1 : 1));
				data[n][i] /= 8388607.0; //2^23-1;
					       
			} else if (BitsPerSample == 32) {
				data[n][i] = ((double) ((int32_t *) block)[m]) / 2147483647.0; // 2^31-1
			} else {
				free(block);
				return -1;
			}

			data[n][i] *= max;
			} else {
				data[n][i] = 0;
				skipped_samples++;
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

// Write to a fffile assumed to be opened by mkff.

int ffwritenext(FFFILE *fffile, FF_data_header *header, double **data, uint32_t n_samples) {

// Again, we assume we're on a little-endian system (like the x86 family of
// processors)

	FILE *file = fffile->fp;
	uint16_t NumChannels = fffile->fmt.NumChannels;
	uint16_t BitsPerSample = fffile->fmt.BitsPerSample;
	void *block = malloc(NumChannels*(BitsPerSample/8));
	uint32_t i,n,m;
	int32_t int24t;
	int32_t skipped_samples;
	int32_t skipped_samples_total = 0;
	uint8_t marker = 42;
	char b = '\0';
	float max;
	long tmp;
	uint32_t size_diff;
	struct stat stat_result;
	double next;
	if (fffile->mode != FFWRITE) {
		fprintf(stderr, "Attempting to write read-only ff file\n");
		return -1;
	}

	if (n_samples%2 != 0) {  // real_part, complex_part, ...
		fprintf(stderr, "Number of samples to fffiles must be even.\n");
		return -1;
	}
	// Assume we're at a new frame at end of file.
	if (fwrite(&marker, 1, 1, file) != 1) return -1;
	if (fwrite(&(header->real_max), sizeof(float), 1, file) != 1) return -1;
	if (fwrite(&(header->imag_max), sizeof(float), 1, file) != 1) return -1;

	for (n = 0; n < NumChannels; n++)
		if (fwrite(header->bitmap[n], header->bitmap_size, 1, file) != 1) return -1;

	for (i = 0; i < n_samples; i++) {
		max = ((i%2 == 0) ? header->real_max : header->imag_max);
		skipped_samples = 0;
		for (n = 0; n < NumChannels; n++) {
		if (bitmap_handler_get(header, i/2, n)) {
		m = n - skipped_samples;
		next = data[n][i]/max;
		if (BitsPerSample == 8) {
			((int8_t *) block)[m] = (int8_t) (next*127); 
		} else if (BitsPerSample == 16) {
			((int16_t *) block)[m] = (int16_t) (next*32767.0);
		} else if (BitsPerSample == 24) {
			int24t = (int32_t) next * 8388607.0;
			((uint8_t *) block)[m] = (uint8_t) (int24t & 0xFF);
			((uint8_t *) block)[m+1] = (uint8_t) (int24t & 0xFF00) >> 8;
			((uint8_t *) block)[m+2] = (uint8_t) (int24t & 0xFF0000) >> 16;
		} else if (BitsPerSample == 32) {
			((int32_t *) block)[m] = (uint32_t) (next*2147483647.0);
		}

		} else {
			skipped_samples++;
		}
		}
		skipped_samples_total += skipped_samples;
		if (skipped_samples != NumChannels)
		if (fwrite(block, (NumChannels-skipped_samples)*(BitsPerSample/8), 1, file) != 1) {
			fprintf(stderr,"Could not write to fffile. It might be corrupt.\n");
			exit(-1);
		}
	}

	free(block);
		// Update sizes
	tmp = ftell(file); // Save old position
	fstat(fileno(file), &stat_result); // Get file size
	if (tmp % 2 == 1) {
		if(fwrite(&b, 1,1,file) != 1) {
			fprintf(stderr,"Could not write header data\n");
			exit(-1);
		}
	}
	size_diff = stat_result.st_size - (fffile->header.ChunkSize + 8);
	fffile->header.ChunkSize += size_diff;
	fffile->data.SubchunkSize += size_diff;
	fseek(file, 4, SEEK_SET); // Where chunkID is stored
	if(fwrite(&(fffile->header.ChunkSize), 4, 1, file) != 1) {
		fprintf(stderr,"Could not write header data\n");
		exit(-1);
	}
	fseek(file, fffile->data.start - 4, SEEK_SET); // Where subchunkID is
	if(fwrite(&(fffile->data.SubchunkSize), 4, 1, file) != 1) {
		fprintf(stderr,"Could not write header data\n");
		exit(-1);
	}
	fseek(file, tmp, SEEK_SET); // Go back to start

	return 0;
}


void bitmap_handler_alloc(FF_data_header *h, uint32_t num_bits, int num_channels) {
	uint32_t s;
	int i;
	s = num_bits/8;
	if (num_bits % 8 != 0) s++;
	
	h->num_channels = num_channels;
	h->bitmap_size = s;
	h->bitmap = (uint8_t **) malloc(num_channels*sizeof(uint8_t *));
	for (i = 0; i < num_channels; i++) {
		h->bitmap[i] = (uint8_t *) malloc(s*sizeof(uint8_t));
	}

}

void bitmap_handler_free(FF_data_header *h) {
	int i;
	for (i = 0; i < h->num_channels; i++)
		free(h->bitmap[i]);
		
	free(h->bitmap);
	h->bitmap_size = 0;
}
void bitmap_handler_set(FF_data_header *h, uint32_t bitno, int channel) {

	(h->bitmap[channel][bitno/8]) |= 1 << (bitno % 8);

}
void bitmap_handler_clr(FF_data_header *h, uint32_t bitno, int channel) {
	(h->bitmap[channel][bitno/8]) &= ~(1 << (bitno % 8));
}

int bitmap_handler_get(FF_data_header *h, uint32_t bitno, int channel) {
	if (h->bitmap[channel][bitno/8] & ((1 << (bitno % 8))))
		return 1;

	return 0;
}

