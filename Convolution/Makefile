ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

# Some of these may be empty; that's fine
ifeq ($(GPU_CONVOLUTION_ENABLED),GPU_CONVOLUTION_ENABLED)
	CXXFLAGS += $(CUDA_LIBS)
endif

CXXFLAGS += $(PARSEHEADER_LIBS) $(THREAD_LIBS) $(LIBS)

# -I $(ROOT_DIR)/singlestep
CPPFLAGS += $(PARSEHEADER_CPPFLAGS) $(THREAD_CPPFLAGS) -I $(ROOT_DIR)/singlestep/Multipoles

ConvolutionDriver: convolutionwrapper.cpp
	$(CXX) -o $@ $< $(CPPFLAGS) $(CXXFLAGS) -MMD

-include ConvolutionDriver.d

distclean: clean
	$(RM) *.d

clean:
	$(RM) *.o ConvolutionDriver

.PHONY: clean distclean
