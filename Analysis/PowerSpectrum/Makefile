CXX = icc
CXXFLAGS = -O3 -fbuiltin -fp-model fast=2  -ip -xHost -mavx -openmp 
CPPFLAGS = -shared -fPIC
INCLUDES = -I../../include/ -I../../Library/include
LIBS = -liomp5 -lgsl -lgslcblas -lfftw3_omp -lfftw3f -lfftw3 -lm -L../../Library/lib -labacus_avx


all: libpowerf.so libpowerd.so

libpowerf.so: power.cc
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)
	
libpowerd.so: power.cc
	$(CXX) ${CXXFLAGS} -DDOUBLE_PRECISION $(CPPFLAGS) -o $@ $^ $(LIBS)


#%.o:%.cpp Makefile
#	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)	

clean:
	rm -f *.o
	
distclean:
	rm -f *.o *.so
	
.PHONY: clean distclean all
