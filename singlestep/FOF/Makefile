CC = g++ -march=native
#CC = icc -xHost
OMP = -DOMP -fopenmp
TBB = -DTBB
CCFLAGS = $(OMP) $(TBB) -I../../include/ -I../Timestep -I../Direct -I../DataModel -I../../ParseHeader -I../Output -I../../Convolution -std=c++11 -O0 -g -Wall
LIBS = -L../../ParseHeader -lparseheader -ltbb

BINARIES = test_driver fof_sublist slab_accum grouplink standalone_fof microstep_test.so

test_driver: test_driver.cpp
	$(CC) $(CCFLAGS) test_driver.cpp $(LIBS) -o test_driver

fof_sublist: fof_sublist.cpp
	$(CC) -DTEST $(CCFLAGS) fof_sublist.cpp -o fof_sublist

slab_accum: slab_accum.cpp
	$(CC) -DTEST $(CCFLAGS) slab_accum.cpp -o slab_accum

grouplink: grouplink.cpp
	$(CC) -DTEST $(CCFLAGS) grouplink.cpp -o grouplink

standalone_fof: dependencies_fof.cpp
	$(CC) $(CCFLAGS) dependencies_fof.cpp $(LIBS) -o standalone_fof

microstep_test: test_driver_microstep.cpp microstep.cpp
	$(CC) $(CCFLAGS) $< -shared -fPIC -o $@.so

default: test_driver

tar:
	tar cvfz groupfinding.tgz *cpp Support/* Makefile readme.groupfinding.txt

.PHONY: clean

clean:
	$(RM) $(BINARIES)
