export CXX = g++ 
export CXXFLAGS = -g -O0 -fopenmp -mavx -Wall \
	-DMAXCPD=8192 -DDOUBLE_PRECISION \
	-I../Direct -I../include -I../Multipoles -I../ParseHeader -I../Convolution

#export FFTW = /Users/eisenste/fftw3 -std=c++0x
export ASMLIB = ../Multipoles/C2R.a ../Multipoles/ETASM.o ../Multipoles/CMASM.o

-include ../Makefile.local

default: singlestep

singlestep: singlestep.o
	$(CXX) $(CXXFLAGS) \
	    -o $@ $< $(ASMLIB)\
	    -L../ParseHeader/ -lparseheader -lfftw3 -lgomp

singlestep.o: singlestep.cpp
	$(CXX) $(CXXFLAGS) -I$(FFTW)/include \
	    -DGITVERSION=\"`git rev-parse HEAD`\" \
	    -c -o $@ $< 

clean:
	rm -f *.o 

distclean:
	rm -f *.o singlestep
