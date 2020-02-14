ROOT_DIR := ../..
-include $(ROOT_DIR)/common.mk

ABACUS_INC := ${ROOT_DIR}/include

CXXFLAGS ?= -O3 -std=c++11
CPPFLAGS ?= -I${ABACUS_INC}

fast_cksum.so: $(ABACUS_INC)/crc32_fast.cpp $(ABACUS_INC)/crc32_fast.h
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -o $@ $< -shared -fPIC
