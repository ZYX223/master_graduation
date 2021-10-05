#include "dmmo.h"

//统计通信量前对内存物理地址空间建立哈希表

static struct dmmo_mem_info *mem;
unsigned long dmmo_addr_conflict;


static inline
unsigned long addr_to_page(unsigned long address)
{
	return (address >> dmmo_shift);
}


static inline
struct dmmo_mem_info* dmmo_get_mem(unsigned long address)
{
	return &mem[hash_32(addr_to_page(address), dmmo_mem_hash_bits)];
}


/* get mem elem, initialize if necessary */
struct dmmo_mem_info* dmmo_get_mem_init(unsigned long address)
{
	struct dmmo_mem_info *elem = dmmo_get_mem(address);
	unsigned long page = addr_to_page(address);

	if (elem->pg_addr != page) { /* new elem */
		if (elem->pg_addr != 0) { /* delete old elem */
			dmmo_addr_conflict++;
		}

		elem->sharer[0] = -1;
		elem->sharer[1] = -1;
		elem->pg_addr = page;
	}

	return elem;
}


void dmmo_mem_init(void)
{
	int hash_size = 1UL << dmmo_mem_hash_bits;

	if (!mem)
		mem = vmalloc(sizeof(struct dmmo_mem_info) * hash_size);

	if (mem)
		memset(mem, 0, sizeof(struct dmmo_mem_info) * hash_size);
	else
		printk("DMMO BUG: could not allocate memory for mem hash table\n");

	dmmo_addr_conflict = 0;
}


void dmmo_mem_cleanup(void)
{
	if (mem)
		vfree(mem);
}
