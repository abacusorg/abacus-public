CC ?= $(shell which icc >/dev/null && echo "icc" || echo "g++")
export CC
all: AbacusCosmo.py

VPATH=../../include

# Try, but don't hard fail because some systems don't have swig
# and this isn't a mandatory module
AbacusCosmo.py: Cosmology.cpp Cosmology.h AbacusCosmo.i setup.py
	-python setup.py build_ext --inplace
    
clean:
	$(RM) -rf build/ AbacusCosmo_wrap.cpp
    
distclean: clean
	$(RM) *.so *.pyc AbacusCosmo.py
	
test:
	pytest

.PHONY: all clean distclean test
