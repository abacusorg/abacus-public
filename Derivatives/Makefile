#CXX = g++ 
CXXFLAGS = -O2 -openmp -mavx -I../include 

CPPFLAGS = -I../include
-include ../Makefile.local

CreateDerivatives: CreateDerivatives.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) \
	-o $@ $< \
	-lfftw3 


clean: 
	rm -f CreateDerivatives *.o