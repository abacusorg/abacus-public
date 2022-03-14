ROOT_DIR := ..
-include $(ROOT_DIR)/common.mk

# allow one to type `make QUAD=1` for convenience
QUAD ?= 0

CPPFLAGS += $(THREAD_CPPFLAGS) -I $(ROOT_DIR)/singlestep/Multipoles -I $(ROOT_DIR)/include -DQUAD_DOUBLE=$(QUAD)
CXXFLAGS += $(THREAD_LIBS) $(SINGLESTEP_LIBS) $(FFTW_OMP_LIBS) -std=gnu++17

CreateDerivatives: CreateDerivatives.cpp
	$(CXX) $(CPPFLAGS) -MMD $< -o $@ $(CXXFLAGS)

-include CreateDerivatives.d

clean:
	$(RM) CreateDerivatives

distclean: clean
	$(RM) *.d

.PHONY: clean distclean
