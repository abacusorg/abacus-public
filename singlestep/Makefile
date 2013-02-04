export CXX = g++ 
export CXXFLAGS = -g -O0 -fopenmp -std=c++0x -Wall \
	-DMAXCPD=8192 \
	-I../Direct -I../include -I../Multipoles -I../ParseHeader -I../Convolution

export FFTW = /Users/eisenste/fftw3
#export ASMLIB = ../Multipoles/C2R.a ../Multipoles/ETASM.o ../Multipoles/CMASM.o

-include ../Makefile.local

default: singlestep

singlestep: singlestep.o
	$(CXX) $(CXXFLAGS) \
	    -o $@ $< \
	    -L../ParseHeader/ -lparseheader -L$(FFTW)/lib -lgomp

singlestep.o: singlestep.cpp
	$(CXX) $(CXXFLAGS) -I$(FFTW)/include \
	    -DGITVERSION=\"`git rev-parse HEAD`\" \
	    -c -o $@ $< 

clean:
	rm -f *.o 

distclean:
	rm -f *.o singlestep
