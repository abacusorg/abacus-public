#export CXX = g++
CXXFLAGS = -mavx -O0 -g3 -openmp -DMAXCPD=81922 -L../ParseHeader $(VERSIONFLAGS)

CPPFLAGS = -I../include -I../Library/include -I../Library/Multipoles -I../singlestep -I../ParseHeader

-include ../Makefile.local

ConvolutionDriver: convolutionwrapper.cpp ../include/Parameters.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) \
	-o $@ $< ../$(ABACUS_VER).a \
	-L$(FFTW)/lib -lfftw3 -lparseheader 

distclean:
	$(RM) -f *.o ConvolutionDriver

clean:
	$(RM) -f *.o

