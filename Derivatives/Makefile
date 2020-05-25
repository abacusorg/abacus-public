ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

CPPFLAGS += $(THREAD_CPPFLAGS) -I $(ROOT_DIR)/singlestep/Multipoles
CXXFLAGS += $(THREAD_LIBS) $(SINGLESTEP_LIBS) $(FFTW_OMP_LIBS)

CreateDerivatives: CreateDerivatives.cpp
	$(CXX) $(CPPFLAGS) -MMD $< -o $@ $(CXXFLAGS)

-include CreateDerivatives.d

clean:
	$(RM) CreateDerivatives

distclean: clean
	$(RM) *.d

.PHONY: clean distclean
