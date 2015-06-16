#CXX = icc -liomp5 -openmp
#CXXFLAGS = -O3 -fbuiltin -fp-model fast=2  -ip -xHost -openmp 
CXX = g++ -lgomp -fopenmp
CXXFLAGS = -O3 -fbuiltin -fopenmp 
CPPFLAGS = -shared -fPIC
INCLUDES = -I../../include/ -I../../Convolution
LIBS = -lgsl -lgslcblas -lfftw3_omp -lfftw3f_omp -lfftw3f -lfftw3 -lm


all: libpowerf.so libpowerd.so

libpowerf.so: power.cpp Makefile
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) $(INCLUDES) -o $@ $< $(LIBS)
	
libpowerd.so: power.cpp Makefile
	$(CXX) ${CXXFLAGS} -DDOUBLEPRECISION $(CPPFLAGS) $(INCLUDES) -o $@ $< $(LIBS)


#%.o:%.cpp Makefile
#	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)	

clean:
	rm -f *.o
	
distclean:
	rm -f *.o *.so
	
.PHONY: clean distclean all
