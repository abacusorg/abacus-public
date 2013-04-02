#CXX = g++ 
CXXFLAGS = -O3 -openmp -mavx -I../include  -I../Library/Multipoles -I../Library/include

CPPFLAGS = -I../include
-include ../Makefile.local

CreateDerivatives: CreateDerivatives.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) \
	-o $@ $< \
	-lfftw3 


clean: 
	rm -f *.o
	
distclean: 
	rm -f CreateDerivatives *.o