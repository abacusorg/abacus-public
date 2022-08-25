ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

all: AbacusCosmoExt

COSMO_INCLUDES = $(addprefix $(ROOT_DIR)/include/, Cosmology.cpp Cosmology.h)

export CC := $(PYTHON_CXX)  # Seems to assume CC can compile C++
export CXX := $(PYTHON_CXX)
export LDFLAGS := $(PYTHON_LDFLAGS)
unexport CPPFLAGS

AbacusCosmoExt: $(COSMO_INCLUDES) setup.py
	python setup.py build_ext --inplace

clean:
	$(RM) -r build/ *.so *.pyc AbacusCosmo.py

distclean: clean
	
test:
	pytest

.PHONY: all clean distclean test AbacusCosmoExt
