CXX = $(shell which icc &>/dev/null && echo "icc" || echo "g++")
CC = $(CXX)
CXXFLAGS = -O1
INCLUDES = -I../include/
LIBS = -lm

CFLAGS = -O1 ${INCLUDES}

TARGET = phheader phdata \
dump_double3 dump_float3 dump_uint64 \
dump_rvdouble dump_pack14 pack14_to_rvdouble

all: ${TARGET}

%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)

clean:
	rm -f *.o
	
distclean:
	rm -f *.o ${TARGET}
	
.PHONY: clean distclean all
