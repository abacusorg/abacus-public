CXX = g++ 
CXXFLAGS = -O2 -fopenmp -mavx -fopenmp -I../include 

FFTW = /Users/eisenste/fftw3

-include ../Makefile.local

CreateDerivatives: CreateDerivatives.cpp
	$(CXX) $(CXXFLAGS) -I$(FFTW)/include \
	-o $@ $< \
	-L$(FFTW)/lib -lfftw3 
