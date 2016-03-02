ifneq ($(DISABLE_GPU),1)

CUDACC=nvcc -O3

ARCHFLAGS=  --generate-code arch=compute_35,code=sm_35\
	    --generate-code arch=compute_30,code=sm_30\
	    --generate-code arch=compute_52,code=sm_52

CUDAFLAGS= $(VERSIONFLAGS) $(ARCHFLAGS) -I$(ABACUS)/include -lineinfo -m64 -Xptxas="-v" -maxrregcount 47


gpudirect.d: 
	$(CUDACC) $(CUDAFLAGS) -E -Xcompiler "-isystem $(CUDA_INSTALL_PATH)/include\
	    -MM" -I$(ABACUS)/include -o $@ GPUPencilDirect.cu

include gpudirect.d

GPUPencilDirect.o: gpudirect.d Makefile ../Makefile
	$(CUDACC) $(CUDAFLAGS) -o $@ -c $(@:.o=.cu)

gpudirect.a: GPUPencilDirect.o
	$(CUDACC) $(CUDAFLAGS) -lib $^ -o $@ 

endif

clean:
	rm -rf *.o *.ptx *.fatbin *.cubin nvtemp

distclean:
	rm -rf *.o *.d *.a *.ptx *.fatbin *.cubin nvtemp DirectBench

