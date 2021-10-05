#ifndef __DMMO_COMM_H
#define __DMMO_COMM_H

#include <linux/delay.h>
#include <linux/hash.h>
#include <linux/kthread.h>
#include <linux/kallsyms.h>
#include <asm-generic/tlb.h>

#include <linux/version.h>

#include "libmapping.h"

#define TSC_DELTA 100000*1000*1000UL

#define DMMO_PID_HASH_BITS 14UL
#define DMMO_PID_HASH_SIZE (1UL << DMMO_PID_HASH_BITS)

struct dmmo_mem_info {
	unsigned long pg_addr;
	unsigned long tsc;
	s16 sharer[2];
};

extern int max_threads;
extern int max_threads_bits;
extern int dmmo_shift;
extern int dmmo_mem_hash_bits;

extern int do_pf;

extern unsigned long dmmo_pte_fixes;
extern unsigned long dmmo_pf;
extern unsigned long dmmo_addr_conflict;
extern unsigned long dmmo_pf_extra;
extern int dmmo_vma_shared_flag;

void reset_stats(void);

/* PID/TID functions */
int dmmo_get_tid(int pid);
int dmmo_get_pid(int tid);
int dmmo_add_pid(int pid);
void dmmo_delete_pid(int pid);
void dmmo_pid_clear(void);
int dmmo_get_num_threads(void);
int dmmo_get_active_threads(void);

/* Mem functions */
struct dmmo_mem_info* dmmo_get_mem_init(unsigned long addr);
void dmmo_mem_init(void);
void dmmo_mem_cleanup(void);

/* Communication functions */
void dmmo_check_comm(int tid, unsigned long address);
void dmmo_print_comm(void);
void dmmo_comm_init(void);

/* Extra pagefaults thread */
int dmmo_pagefault_func(void* v);
void dmmo_pf_thread_init(void);
extern int dmmo_num_faults;

/* Mapping */
int dmmo_map_func(void* v);
void dmmo_set_affinity(int tid, int core);
void dmmo_map_init(void);

/* Probes */
void dmmo_probes_init(void);
void dmmo_probes_cleanup(void);

/* Intercept */
int interceptor_start(void);
void interceptor_stop(void);

/* Topology */
void topo_init(void);
extern int num_nodes, num_cpus, num_cores, num_threads, pu[256];

/* Communication Matrix */
extern struct dmmo_comm_matrix dmmo_matrix;

static inline
unsigned get_comm(int first, int second)
{
	if (first > second)
		return dmmo_matrix.matrix[(first << max_threads_bits) + second];
	else
		return dmmo_matrix.matrix[(second << max_threads_bits) + first];
}


static inline
void inc_comm(int first, int second, unsigned old_tsc, unsigned long new_tsc)
{
	/* TODO:
		- replace with Atomic incr.
		- check if spinlock needed
		- compare TSC_DELTA
	*/

	if (first > second)
		dmmo_matrix.matrix[(first << max_threads_bits) + second] ++;
	else
		dmmo_matrix.matrix[(second << max_threads_bits) + first] ++;
}

/* ProcFS */
int dmmo_proc_init (void);
void dmmo_proc_cleanup(void);


#endif
