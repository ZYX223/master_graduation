#include "dmmo.h"
#include "version.h"

#include <linux/module.h>

MODULE_LICENSE("GPL");

#define NUM_FAULTS_DEFAULT 9
#define NUM_MAX_THREADS_DEFAULT 1024
#define DMMO_SHIFT_DEFAULT 12
#define DMMO_MEM_HASH_BITS_DEFAULT 22

#ifdef ENABLE_EXTRA_PF
	#define ENABLE_EXTRA_PF 1
#else
	#define ENABLE_EXTRA_PF 0
#endif

static struct task_struct *dmmo_pf_thread;
static struct task_struct *dmmo_map_thread;

//声明属性参数
int num_faults = NUM_FAULTS_DEFAULT;
int do_map = 0;
int max_threads = NUM_MAX_THREADS_DEFAULT;
int max_threads_bits = 0;
int dmmo_shift = DMMO_SHIFT_DEFAULT;
int dmmo_mem_hash_bits = DMMO_MEM_HASH_BITS_DEFAULT;
int do_pf = ENABLE_EXTRA_PF;        

//编写内核模块时通过module_param()向当前模块传递属性参数
module_param(num_faults, int, 0);
module_param(do_map, int, 0);
module_param(max_threads, int, 0);
module_param(dmmo_shift, int, 0);
module_param(dmmo_mem_hash_bits, int, 0);
module_param(do_pf, int, 0);

//初始化统计矩阵
struct dmmo_comm_matrix dmmo_matrix = {.matrix = NULL, .nthreads = 0};

int dmmo_num_faults;

//初始化模块
int init_module(void)
{
	max_threads = roundup_pow_of_two(max_threads);   
	max_threads_bits = ilog2(max_threads);    
	dmmo_num_faults = num_faults;   

	//printk("DMMO: Start (version %s)\n", DMMO_VERSION);
	//打印映射机制中的一些属性信息
	printk("DMMO: Start\n");
	printk("    additional pagefaults (do_pf): %s %s\n", do_pf ? "yes" : "no", do_pf==ENABLE_EXTRA_PF ? "(default)" : "");
	printk("    extra pagefaults (num_faults): %d %s\n", num_faults, num_faults==NUM_FAULTS_DEFAULT ? "(default)" : "");
	printk("    maximum threads (max_threads): %d %s\n", max_threads, max_threads==NUM_MAX_THREADS_DEFAULT ? "(default)" : "");
	printk("    use mapping (do_map): %s\n", do_map ? "yes" : "no (default)");
	printk("    granularity (dmmo_shift): %d bits %s\n", dmmo_shift, dmmo_shift==DMMO_SHIFT_DEFAULT ? "(default)" : "");
	printk("    mem hash table size (dmmo_mem_hash_bits): %d bits, %d elements %s\n", dmmo_mem_hash_bits, 1<<dmmo_mem_hash_bits, dmmo_mem_hash_bits==DMMO_MEM_HASH_BITS_DEFAULT ?
		"(default)" : "");

	spin_lock_init(&dmmo_matrix.lock);//保护矩阵数据

	reset_stats();//重置属性设置
	dmmo_probes_init();//内核探针

	dmmo_proc_init();//绑定设置初始化
	dmmo_map_init();//映射结果初始化

    //如果引入extra page fault，创建一个线程处理额外页错误
	if (do_pf) {
		dmmo_pf_thread = kthread_create(dmmo_pagefault_func, NULL, "dmmo_pf_thread");  //kthread_create创建线程后，不会马上运行，而是将该函数返回的task_struct类型指针
																                       //dmmo_pf_thread传递给wake_up_process()，然后通过此函数运行线程。
		wake_up_process(dmmo_pf_thread);
	}

	topo_init();

	//如果执行映射，创建一个线程处理映射
	if (do_map) {
		dmmo_map_thread = kthread_create(dmmo_map_func, NULL, "dmmo_map_thread");
		wake_up_process(dmmo_map_thread);
	}

	return 0;
}

//清除模块
void cleanup_module(void)
{
	if (dmmo_pf_thread)
		kthread_stop(dmmo_pf_thread);

	if (dmmo_map_thread)
		kthread_stop(dmmo_map_thread);

	if (dmmo_matrix.matrix)
		kfree(dmmo_matrix.matrix);

	dmmo_proc_cleanup();

	dmmo_probes_cleanup();

	dmmo_mem_cleanup();

	//printk("DMMO: Quit (version %s)\n", DMMO_VERSION);
	printk("DMMO: Quit\n");
}


