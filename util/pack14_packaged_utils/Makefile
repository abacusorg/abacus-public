CPPFLAGS ?= -I include/
CPPFLAGS += -O3

TARGET := phheader phdata \
dump_pack14 pack14_to_rvdouble pack14_to_rvfloat

SAMPLE_FILE := sample.pack14.dat

all: ${TARGET}

%.o:%.cpp
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)

test: dump_pack14 phdata
	./phdata $(SAMPLE_FILE) | ./dump_pack14 -

clean:
	$(RM) ${TARGET} *.o
	
distclean: clean
	
.PHONY: clean distclean all test
