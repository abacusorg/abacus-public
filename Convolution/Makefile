ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

ifeq ($(GPU_CONVOLUTION_ENABLED),GPU_CONVOLUTION_ENABLED)
	CXXFLAGS += $(CUDA_LIBS)
endif

CXXFLAGS += $(PARSEHEADER_LIBS) $(THREAD_LIBS) $(SINGLESTEP_LIBS)

CPPFLAGS += $(PARSEHEADER_CPPFLAGS) $(THREAD_CPPFLAGS) $(TCMALLOC_CPPFLAGS) -I $(ROOT_DIR)/singlestep/Multipoles -I $(ROOT_DIR)/DirectIO

convolution: convolutionwrapper.cpp
	$(CXX) $(CPPFLAGS) -o $@ $< -MMD $(CXXFLAGS)

-include convolution.d

distclean: clean
	$(RM) *.d

clean:
	$(RM) *.o convolution

.PHONY: clean distclean
