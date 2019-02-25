ROOT_DIR := ../..

# We invoke the "detect_gpus" recipe from the ./configure script, but common.mk won't yet exist!
ifneq ($(MAKE_CONTEXT),AUTOCONF)
include $(ROOT_DIR)/common.mk
endif

ifeq ($(CUDA_ENABLED),CUDA_ENABLED)

ARCHFLAGS = --generate-code arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH)\

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

# Don't want make distclean to force dependency generation
# Only have to do this because nvcc doesn't offer -MMD to do inline dep generation
ifeq (,$(findstring clean,$(MAKECMDGOALS)))
-include GPUControl.d
endif

GPUControl.d:
	$(NVCC) $(NVCC_FLAGS) -Wno-deprecated-gpu-targets -Xcompiler "$(NVCC_XCOMPILERFLAGS)" -M $(@:.d=.cu) -o $@

GPUControl.o: GPUControl.cu GPUControl.d
	$(NVCC) $(CUDAFLAGS) -Xcompiler "$(NVCC_XCOMPILERFLAGS)" -o $@ -c $<

else
all:
endif

BINARIES=avx512direct detect_gpus

# NVCC may not exist outside the CUDA_ENABLED region!
# But ./configure won't invoke this unless we do have CUDA
detect_gpus: detect_gpus.cu
	$(NVCC) -Wno-deprecated-gpu-targets -o $@ $< -ccbin $(CXX)

avx512direct: avx512direct.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DTEST -o $@ $<

clean:
	$(RM) *.o *.ptx *.fatbin *.cubin nvtemp $(BINARIES) *.a

distclean: clean
	$(RM) *.d

.PHONY: clean distclean all
