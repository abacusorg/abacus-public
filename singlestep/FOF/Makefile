CC = g++
OMP = -DOMP -fopenmp -lgomp 
TBB = -DTBB
CCFLAGS = $(OMP) $(TBB) -I../../include/ -I../Timestep -I../Direct -std=c++11 -O3 -g -march=native

test_driver: test_driver.cpp
	$(CC) $(CCFLAGS) test_driver.cpp -ltbb -o test_driver

fof_sublist: fof_sublist.cpp
	$(CC) -DTEST $(CCFLAGS) fof_sublist.cpp -o fof_sublist

slab_accum: slab_accum.cpp
	$(CC) -DTEST $(CCFLAGS) slab_accum.cpp -o slab_accum

grouplink: grouplink.cpp
	$(CC) -DTEST $(CCFLAGS) grouplink.cpp -o grouplink

dependencies_fof: dependencies_fof.cpp
	$(CC) $(CCFLAGS) dependencies_fof.cpp -ltbb -o dependencies_fof

microstep_test: test_driver_microstep.cpp microstep.cpp
	$(CC) $(CCFLAGS) $< -shared -fPIC -o $@.so
	#python build_microstep_ffi.py

default: test_driver

tar:
	tar cvfz groupfinding.tgz *cpp Support/* Makefile readme.groupfinding.txt
