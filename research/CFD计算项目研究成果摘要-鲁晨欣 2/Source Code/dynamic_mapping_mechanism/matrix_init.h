#ifndef LIBDMMO_H
#define LIBDMMO_H

#define DMMO_MAGIC 0x48151623

#ifdef __KERNEL__
	#include <linux/slab.h>
#else
	#include <string.h>
	#include <stdlib.h>
	#include <stdio.h>
#endif

typedef struct dmmo_matrix{
	int num_threads;
	int max_threads;
	int *pids;
	unsigned *matrix;
} dmmo_matrix_t;

void *dmmo_matrix_encode(dmmo_matrix_t *data);
dmmo_matrix_t dmmo_matrix_decode(void *ptr);
size_t dmmo_matrix_size(void *ptr);

#ifndef __KERNEL__
unsigned *dmmo_get_small_matrix(dmmo_matrix_t *data);
void dmmo_matrix_print(dmmo_matrix_t *data);
#endif

#endif /* end of include guard: LIBDMMO_H */
