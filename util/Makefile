#CXX = icc
CXXFLAGS = -O1
INCLUDES = -I../../include/ -lm
LIBS = 

CFLAGS = -O1 ${INCLUDES}

TARGET = phheader phdata \
dump_double3 dump_float3 dump_uint64 \
dump_rvdouble dump_pack14

all: ${TARGET}

%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)

clean:
	rm -f *.o
	
distclean:
	rm -f *.o ${TARGET}
	
.PHONY: clean distclean all
