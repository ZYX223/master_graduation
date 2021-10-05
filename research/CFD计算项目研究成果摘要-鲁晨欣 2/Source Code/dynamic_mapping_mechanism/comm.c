#include "dmmo.h"

//通信量矩阵更新

static inline
int get_num_sharers(struct dmmo_mem_info *elem)
{
	if (elem->sharer[0] == -1 && elem->sharer[1] == -1)
		return 0;
	if (elem->sharer[0] != -1 && elem->sharer[1] != -1)
		return 2;
	return 1;
}


void dmmo_check_comm(int tid, unsigned long address)
{
	struct dmmo_mem_info *elem = dmmo_get_mem_init(address);
	int sh = get_num_sharers(elem);
	unsigned new_tsc = get_cycles();

	switch (sh) {
		case 0: /* no one accessed page before, store accessing thread in pos 0 */
			elem->sharer[0] = tid;
			break;

		case 1: /* one previous access => needs to be in pos 0 */
			if (elem->sharer[0] != tid) {
				inc_comm(tid, elem->sharer[0], elem->tsc, new_tsc);
				elem->sharer[1] = elem->sharer[0];
				elem->sharer[0] = tid;
			}
			break;

		case 2: /* two previous accesses */
			if (elem->sharer[0] != tid && elem->sharer[1] != tid) {
				inc_comm(tid, elem->sharer[0], elem->tsc, new_tsc);
				inc_comm(tid, elem->sharer[1], elem->tsc, new_tsc);
				elem->sharer[1] = elem->sharer[0];
				elem->sharer[0] = tid;
			} else if (elem->sharer[0] == tid) {
				inc_comm(tid, elem->sharer[1], elem->tsc, new_tsc);
			} else if (elem->sharer[1] == tid) {
				inc_comm(tid, elem->sharer[0], elem->tsc, new_tsc);
				elem->sharer[1] = elem->sharer[0];
				elem->sharer[0] = tid;
			}

			break;
	}

	elem->tsc = new_tsc;
}


void dmmo_print_comm(void)
{
	int i, j;
	int nt = dmmo_get_num_threads();
	unsigned long sum = 0, sum_sqr = 0;
	unsigned long av, va;

	if (nt < 2)
		return;

	for (i = nt-1; i >= 0; i--) {
		for (j = 0; j < nt; j++) {
			int s = get_comm(i,j);
			sum += s;
			sum_sqr += s*s;
			printk ("%u", s);
			if (j != nt-1)
				printk (",");
		}
		printk("\n");
	}

	av = sum / nt / nt;
	va = (sum_sqr - ((sum*sum)/nt/nt))/(nt-1)/(nt-1);

	printk ("avg: %lu, var: %lu, hf: %lu\n", av, va, av ? va/av : 0);
}


void dmmo_comm_init(void)
{
	if (!dmmo_matrix.matrix){
		dmmo_matrix.matrix =  kmalloc (sizeof(unsigned) * max_threads * max_threads, GFP_KERNEL);
		dmmo_matrix.nthreads = 0;
	}

	if (dmmo_matrix.matrix)
		memset(dmmo_matrix.matrix, 0, sizeof(unsigned) * max_threads * max_threads);
	else
		printk("DMMO BUG: could not allocate memory for comm matrix\n");
}
