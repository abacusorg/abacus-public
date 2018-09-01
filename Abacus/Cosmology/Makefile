ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

all: AbacusCosmo.py

COSMO_INCLUDES = $(addprefix $(ROOT_DIR)/include/, Cosmology.cpp Cosmology.h)

# Export the C++ compiler to swig, which uses the CC variable
# May need to monkey around with other flags if swig flags are wrong for non-GNU compilers
export CC = $(CXX)

# Try, but don't hard fail because some systems don't have swig
# and this isn't a mandatory module
AbacusCosmo.py: $(COSMO_INCLUDES) AbacusCosmo.i setup.py
	python setup.py build_ext --inplace

clean:
	$(RM) -r build/ AbacusCosmo_wrap.cpp *.so *.pyc AbacusCosmo.py

distclean: clean
	
test:
	pytest

.PHONY: all clean distclean test
