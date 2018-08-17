ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

ifeq ($(CUDA_ENABLED),CUDA_ENABLED)

# Weird that CUDA makes you do this manually!
ARCHFLAGS = --generate-code arch=compute_35,code=sm_35\
	    --generate-code arch=compute_30,code=sm_30\
	    --generate-code arch=compute_52,code=sm_52\
	    --generate-code arch=compute_52,code=sm_61\
	    --generate-code arch=compute_60,code=sm_60

# CUDA doesn't like to support recent g++ versions because they're untested,
# but in practice they often work just fine.  The following is an example of
# how to lie to CUDA about the g++ version.
# If CUDA <= 9 and GCC > 6, then override CUDA's protests
# WARNING: for these specific versions, CUDA appears to be truly incompatible,
# hence we will leave it commented
#CUDA_LE_9 := $(shell [ $(CUDA_MAJOR_VER) -le 9 ] && echo true)
#ifeq ($(CUDA_LE_9), true)
#	ifeq ($(CXX),g++)
#		GNUC_GT_6 := $(shell [ $(CXX_MAJOR_VER) -gt 6 ] && echo true)
#		ifeq ($(GNUC_GT_6),true)
#			CPPFLAGS += -D__GNUC__=6
#		endif
#	endif
#endif

# _FORCE_INLINES fixes a weird compilation bug
CPPFLAGS += -D_FORCE_INLINES

CUDAFLAGS = $(ARCHFLAGS) -lineinfo -Xptxas="-v" -maxrregcount 63 --default-stream per-thread $(NVCC_FLAGS)

all: gpudirect.a

gpudirect.a: GPUControl.o
	$(NVCC) $(CUDAFLAGS) -Xcompiler "$(NVCC_XCOMPILERFLAGS)" -lib $^ -o $@ 

-include gpudirect.d

GPUControl.o: GPUControl.cu
	$(NVCC) $(CUDAFLAGS) -Xcompiler "$(NVCC_XCOMPILERFLAGS) -MMD" -o $@ -c $<

else
all:
endif

BINARIES=avx512direct

avx512direct: avx512direct.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DTEST -o $@ $<

clean:
	$(RM) *.o *.ptx *.fatbin *.cubin nvtemp $(BINARIES) *.a

distclean: clean
	$(RM) *.d

.PHONY: clean distclean all
