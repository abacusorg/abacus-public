CXX = g++
CXXFLAGS = -O0 -g3 -fopenmp -DDOUBLE_PRECISION
CPPFLAGS = 
INCLUDES = -I../../include/
LIBS = -lfftw3

all:makespiralics makeanalytic

makespiralics: spiralic.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)

makeanalytic: AnalyticSpiral.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)

%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)
	
	
clean:
	rm -f *.o

distclean:
	rm -f *.o makespiralics makeanalytic analytic
	
.PHONY: clean distclean all
