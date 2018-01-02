CC = g++
OMP = -DOMP -fopenmp -lgomp 
CCFLAGS = $(OMP) -I../../include/ -I../Timestep -std=c++11 -g -mavx

test_driver: test_driver.cpp
	$(CC) $(CCFLAGS) test_driver.cpp -o test_driver

fof_sublist: fof_sublist.cpp
	$(CC) -DTEST $(CCFLAGS) fof_sublist.cpp -o fof_sublist

slab_accum: slab_accum.cpp
	$(CC) -DTEST $(CCFLAGS) slab_accum.cpp -o slab_accum

grouplink: grouplink.cpp
	$(CC) -DTEST $(CCFLAGS) grouplink.cpp -o grouplink

dependencies_fof: dependencies_fof.cpp
	$(CC) $(CCFLAGS) dependencies_fof.cpp -o dependencies_fof

default: test_driver

tar:
	tar cvfz groupfinding.tgz *cpp Support/* Makefile readme.groupfinding.txt
