ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

#CC = icc -xHost
OMP :=  -fopenmp
CPPFLAGS += $(PARSEHEADER_CPPFLAGS) $(OMP) -DOMP -DTBB -I../Timestep -I../Direct -I../DataModel -I../Output -I../../Convolution -std=c++17 -O0
LIBS := $(PARSEHEADER_LIBS)

BINARIES = test_driver fof_sublist slab_accum grouplink standalone_fof

test_driver: test_driver.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test_driver.cpp $(LIBS) -o test_driver

fof_sublist: fof_sublist.cpp
	$(CXX) $(CXXFLAGS) -DTEST $(CPPFLAGS) fof_sublist.cpp $(LIBS) -o fof_sublist

slab_accum: slab_accum.cpp
	$(CXX) $(CXXFLAGS) -DTEST $(CPPFLAGS) slab_accum.cpp $(LIBS) -o slab_accum

grouplink: grouplink.cpp
	$(CXX) $(CXXFLAGS) -DTEST $(CPPFLAGS) grouplink.cpp $(LIBS) -o grouplink

# Note: this interface is deprecated, use abacus/Abacus/standalone_fof.py
standalone_fof: dependencies_fof.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dependencies_fof.cpp $(LIBS) -o standalone_fof

default: test_driver

tar:
	tar cvfz groupfinding.tgz *cpp Support/* Makefile readme.groupfinding.txt

.PHONY: clean default tar

clean:
	$(RM) $(BINARIES) *.o *.so
