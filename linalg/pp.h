#ifndef _PP
#define _PP
#include <pthread.h>
#include <stdlib.h>

// Returns number of logical cpu's on the system.
size_t numCPUs();

typedef struct {
	pthread_mutex_t count_lock;
	pthread_cond_t ok_to_proceed;
	size_t count;
	size_t n_threads;
} barrier_t;

void barrier_init(barrier_t *b, size_t n_threads);
void barrier_destroy(barrier_t *b);
void barrier(barrier_t *b);

#endif
