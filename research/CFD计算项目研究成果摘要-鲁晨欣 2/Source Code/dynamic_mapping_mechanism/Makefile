SRCDIR:=$(dir $(lastword $(MAKEFILE_LIST)))

obj-m := dmmo.o
dmmo-objs := matrix_init.o pagefault.o mem.o pid.o comm.o probes.o map.o topo.o lhmapping.o dmmo_main.o

ccflags-y += -g -Wall -D_DMMO -I$(SRCDIR)/lhmapping $(LMCFLAGS)

isnuma=$(shell hwloc-info|grep -qi numanode; if [ $$? -eq 1 ]; then echo 0; else echo 1;fi)
kmaj=$(shell uname -r | tr '.' ' ' | awk '{print $$1}' )
kmin=$(shell uname -r | tr '.' ' ' | awk '{print $$2}' )
islow=$(shell if [ $(kmaj) -lt 3 -o $(kmin) -lt 8 ]; then echo 1; else echo 0; fi)
addpf=$(shell if [ $(isnuma) -eq 0 -o $(islow) -eq 1 ]; then echo 1; else echo 0; fi)

ifeq ($(addpf), 1)
       ccflags-y += -DENABLE_EXTRA_PF
endif

.PHONY: all clean install dist

all:
	@sync
	@echo "#define DMMO_VERSION \"$(shell git describe); $(shell date)\"" > version.h
	@if stat -t obj/* >/dev/null 2>&1; then mv -f obj/* obj/.*.cmd . ; else mkdir -p obj; fi
	make -j -C /lib/modules/$(shell uname -r)/build M=$(PWD) modules
	@mv -f *.o modules.order Module.symvers dmmo.mod.c .*.cmd obj

clean:
	@rm -rf obj/
	@rm -f version.h
	@rm -f dmmo-*.tar.gz
	make -C /lib/modules/$(shell uname -r)/build M=$(PWD) clean


install: all
	@if test $(shell lsmod | grep dmmo >/dev/null; echo $$?) -eq 0; then sudo /sbin/rmmod dmmo; fi
	sudo /sbin/insmod $(PWD)/dmmo.ko $(options)
	@dmesg | grep -i "dmmo bug"; echo -n

dist: clean
	tar czf dmmo-$(shell git describe).tar.gz *
