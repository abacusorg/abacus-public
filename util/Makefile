#CXX = icc
CXXFLAGS = -O1
INCLUDES = -I../../include/ -lm
LIBS = 

CFLAGS = -O1 ${INCLUDES}

TARGET = phheader phdata dump_double3 dump_rvdouble

all: ${TARGET}

%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)

clean:
	rm -f *.o
	
distclean:
	rm -f *.o ${TARGET}
	
.PHONY: clean distclean all
