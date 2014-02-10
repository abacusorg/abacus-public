export CXX = icc
CXXFLAGS = -mavx -O3 -openmp -DMAXCPD=8192 -L../ParseHeader $(VERSIONFLAGS)

CPPFLAGS = -I../include -I../Library/include -I../Library/lib/Multipoles -I../singlestep -I../ParseHeader

-include ../Makefile.local

ConvolutionDriver: convolutionwrapper.cpp ../include/Parameters.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) \
	-o $@ $< -L../Library/lib -l$(ABACUS_VER) \
	-L$(FFTW)/lib -lfftw3 -lparseheader -liomp5

distclean:
	$(RM) -f *.o ConvolutionDriver

clean:
	$(RM) -f *.o

