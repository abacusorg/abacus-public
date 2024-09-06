ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

ifeq ($(CUDA_ENABLED),CUDA_ENABLED)

ARCHFLAGS := -arch=native

# _FORCE_INLINES fixes a weird compilation bug
CPPFLAGS += -D_FORCE_INLINES

CUDAFLAGS := $(ARCHFLAGS) -lineinfo -Xptxas="-v" --default-stream per-thread $(NVCC_FLAGS)

# for fmt "warning #128-D: loop is not reachable"
CUDAFLAGS += -diag-suppress 128

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

BINARIES=avx512direct

avx512direct: avx512direct.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DTEST -o $@ $<

clean:
	$(RM) *.o *.ptx *.fatbin *.cubin nvtemp $(BINARIES) *.a

distclean: clean
	$(RM) *.d

.PHONY: clean distclean all
