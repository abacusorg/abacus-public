export CXX = g++ 
export CXXFLAGS = -g -O0 -fopenmp  -std=c++0x -Wall  -I../avx -I../include -I../Multipoles -I../singlestep  -DMAXCPD=8192 -I../ParseHeader
#export FFTWLIB = -L/home/drdre/gccfftw3/lib 
#export ASMLIB = ../Multipoles/C2R.a ../Multipoles/ETASM.o ../Multipoles/CMASM.o

-include ../Makefile.local

singlestep: singlestep.o
	$(CXX) $(CXXFLAGS) -o $@ $< -L../ParseHeader/ -lparseheader -lgomp
	

clean:
	rm -f *.o 

distclean:
	rm -f *.o singlestep

default: singlestep
