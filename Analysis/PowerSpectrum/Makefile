ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

CPPFLAGS += -I $(ROOT_DIR)/DirectIO $(THREAD_CPPFLAGS) $(READABACUS_CPPFLAGS) $(GSL_CPPFLAGS)
LDFLAGS = $(CXXFLAGS) -shared $(THREAD_LIBS) $(GSL_LIBS)
CXXFLAGS += -fPIC

unexport CPPFLAGS

# We can't use SINGLESTEP_LIBS; may contain "unsafe" libraries like tcmalloc (see common.mk)
LIBS = $(FFTW_OMP_LIBS) $(READABACUS_LIBS)

PSLIB_SO = PowerSpectrum/pslib$(shell python-config --extension-suffix)

all: lib/libpowerd.so lib/libpowerf.so $(PSLIB_SO)

# TODO: pass the detected FFTW libs from configure to here
lib/libpowerf.so: PowerSpectrum/power.cpp
	mkdir -p lib
	$(CXX) ${LDFLAGS} $(CPPFLAGS) -o $@ $< $(LIBS) -lfftw3f_omp
	
lib/libpowerd.so: PowerSpectrum/power.cpp
	mkdir -p lib
	$(CXX) ${LDFLAGS} -DDOUBLEPRECISION $(CPPFLAGS) $(LIBFLAGS) $(INCLUDES) -o $@ $< $(LIBS) -lfftw3_omp

$(PSLIB_SO):
	CC="$(PYTHON_CXX)" CXX="$(PYTHON_CXX)" LDFLAGS="$(PYTHON_LDFLAGS)" CPPFLAGS="$(PYTHON_CPPFLAGS)" python setup.py build_ext --inplace

clean:
	$(RM) PowerSpectrum/*.o PowerSpectrum/*.so lib/*
	
distclean: clean
	
.PHONY: clean distclean all
