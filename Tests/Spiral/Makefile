CXX = icc
CXXFLAGS = -O0 -g3 -openmp -DDOUBLE_PRECISION
CPPFLAGS = 
INCLUDES = -I../../include/
LIBS = -liomp5

all:makespiralics

makespiralics: spiralic.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)

%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)
	
distclean:
	rm -f *.o makespiralics
	
.PHONY: clean distclean all
