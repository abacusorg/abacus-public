export CXX = g++
export CXXFLAGS = -mavx -O3 -fopenmp -I../include -I../Multipoles -DMAXCPD=81922 -fopenmp

FFTW = /Users/eisenste/fftw3

-include ../Makefile.local

ConvolutionDriver: ConvolutionDriver.cpp
	$(CXX) $(CXXFLAGS) -I$(FFTW)/include \
	-o $@ $< \
	-L$(FFTW)/lib -lfftw3 

distclean:
	$(RM) -f *.o ConvolutionDriver

clean:
	$(RM) -f *.o

