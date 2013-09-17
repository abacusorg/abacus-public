CXX = icc
CXXFLAGS = -O3  -fbuiltin -fp-model fast=2  -ip -xHost -mavx -openmp 
CPPFLAGS = -shared -fPIC
INCLUDES = -I../../include/ -I../../Library/include -I../../Library/lib/convolution -I../../Library/lib/common
LIBS = -liomp5 -lgsl -lgslcblas -lfftw3_omp -lfftw3f_omp -lfftw3f -lfftw3 -lm


all: libpowerf.so libpowerd.so

libpowerf.so: power.cpp Makefile
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) $(INCLUDES) -o $@ $< $(LIBS)
	
libpowerd.so: power.cpp Makefile
	$(CXX) ${CXXFLAGS} -DDOUBLE_PRECISION $(CPPFLAGS) $(INCLUDES) -o $@ $< $(LIBS)


#%.o:%.cpp Makefile
#	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)	

clean:
	rm -f *.o
	
distclean:
	rm -f *.o *.so
	
.PHONY: clean distclean all
