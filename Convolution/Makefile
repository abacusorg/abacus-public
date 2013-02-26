#export CXX = g++
#export CXXFLAGS = -mavx -O3 -fopenmp -DMAXCPD=81922 -fopenmp

CPPFLAGS = -I../include -I../Multipoles

-include ../Makefile.local

ConvolutionDriver: ConvolutionDriver.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) \
	-o $@ $< \
	-L$(FFTW)/lib -lfftw3 

distclean:
	$(RM) -f *.o ConvolutionDriver

clean:
	$(RM) -f *.o

