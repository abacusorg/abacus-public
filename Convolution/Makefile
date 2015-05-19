#export CXX = icc -openmp -liomp5
export CXX = g++ -fopenmp -lgomp
CUDAFLAGS = -DGPUFFT  -lcufft -lcudart -L/usr/local/cuda/lib64 -I/usr/local/cuda/include
FFTWFLAGS = -lfftw3 -L$(FFTW)/lib 

VERSIONFLAGS = -mavx -Ofast $(CUDAFLAGS)
CXXFLAGS = -DMAXCPD=8192 -L../ParseHeader -lparseheader $(VERSIONFLAGS)

CPPFLAGS = -MMD -I../include -I../singlestep/Multipoles -I../singlestep -I../ParseHeader

-include ../Makefile.local

ConvolutionDriver: convolutionwrapper.cpp Makefile
	$(CXX) -o $@ $< $(CXXFLAGS) $(CPPFLAGS)


-include ConvolutionDriver.d

distclean:
	$(RM) -f *.o *.d ConvolutionDriver

clean:
	$(RM) -f *.o *.d

