#export CXX = icc -openmp -liomp5
export CXX = g++ -fopenmp -lgomp
ifneq ($(DISABLE_GPU),1)
CUDAFLAGS = #-DGPUFFT  -lcufft -lcudart -L/usr/local/cuda/lib64 -I/usr/local/cuda/include
endif
FFTWFLAGS = -lfftw3 -lm -L$(FFTW)/lib 

VERSIONFLAGS = -march=native -mavx -O3 $(CUDAFLAGS) $(FFTWFLAGS)
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

