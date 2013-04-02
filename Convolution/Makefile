#export CXX = g++
CXXFLAGS = -mavx -O3 -openmp -DMAXCPD=81922

CPPFLAGS = -I../include -I../Library/include -I../Library/Multipoles -I../singlestep

-include ../Makefile.local

ConvolutionDriver: convolutionwrapper.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) \
	-o $@ $< \
	-L$(FFTW)/lib -lfftw3 

distclean:
	$(RM) -f *.o ConvolutionDriver

clean:
	$(RM) -f *.o

