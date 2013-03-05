CXX = icc
CXXFLAGS = -O3 -openmp -DDOUBLE_PRECISION
CPPFLAGS = -shared -fPIC
INCLUDES = -I../../include/
LIBS = -liomp5 -lgsl -lgslcblas -lm -lfftw3


all: libpower.so

libpower.so: power.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)


%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)	

clean:
	rm -f *.o
	
distclean:
	rm -f *.o *.so
	
.PHONY: clean distclean all