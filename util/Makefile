# touch a FREEZE file to prevent the Makefile from doing anything
# Useful for mixed architecture platforms like Summit and Rhea
FREEZE := $(wildcard FREEZE)

ifneq (,$(FREEZE))
all:
else

ROOT_DIR := ..
-include $(ROOT_DIR)/common.mk

CXX ?= g++
CXXFLAGS ?= -Wall 

CPPFLAGS ?= -I$(ROOT_DIR)/include/
CPPFLAGS += -O3
LIBS ?= -lm

TARGET = phheader phdata \
dump_double3 dump_float3 dump_uint64 \
dump_rvdouble dump_pack14 pack14_to_rvdouble pack14_to_rvfloat

all: ${TARGET}

%.o:%.cpp
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)

clean:
	$(RM) *.o
	
distclean: clean
	$(RM) ${TARGET}
	
.PHONY: clean distclean all

endif  #freeze