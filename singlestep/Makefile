export CXX = g++ 
export CXXFLAGS = -O2  -Wall -fopenmp  -I../avx -I../include -I../Multipoles -I../singlestep  -DMAXCPD=8192 -I../ParseHeader
#export FFTWLIB = -L/home/drdre/gccfftw3/lib 
#export ASMLIB = ../Multipoles/C2R.a ../Multipoles/ETASM.o ../Multipoles/CMASM.o

singlestep: singlestep.o
	$(CXX) $(CXXFLAGS) -o $@ $< -L../ParseHeader/ -lparseheader

clean:
	rm -f *.o singlestep

default: singlestep
