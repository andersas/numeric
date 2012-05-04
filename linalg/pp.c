#include <unistd.h>
#include <pthread.h>
#include "pp.h"


size_t numCPUs() {
 return (size_t) sysconf( _SC_NPROCESSORS_ONLN );
}

// Barriers wait for all threads to reach a certain point
// of evaluation

void barrier_init(barrier_t *b, size_t n_threads) {
	pthread_mutex_init(&(b->count_lock), NULL);
	pthread_cond_init(&(b->ok_to_proceed), NULL);
	b->count = 0;
	b->n_threads = n_threads;
}

void barrier_destroy(barrier_t *b) {
	pthread_mutex_destroy(&(b->count_lock));
	pthread_cond_destroy(&(b->ok_to_proceed));
}

void barrier(barrier_t *b) {

	pthread_mutex_lock(&(b->count_lock));
	b->count++;
	if (b->count == b->n_threads) {
		b->count = 0;
		pthread_cond_broadcast(&(b->ok_to_proceed));
	} else {
		// Spurious wakeups may happen. Reevaluate predicate
		// before sleeping again.
		while (pthread_cond_wait(&(b->ok_to_proceed), &(b->count_lock)) != 0);
	}
	pthread_mutex_unlock(&(b->count_lock));
}

