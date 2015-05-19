CUDACC=nvcc
CUDAFLAGS= -I/usr/local/cuda/include -m64  -arch compute_35\
	   -code sm_35 -lineinfo  -maxrregcount=64 -Xptxas="-v" -O3\
	   -DDIRECTSPLINE
CUDASRC= gpudirect.cu

gpudirect: gpudirect.o gpudirect.d

gpudirect.d:
	$(CUDACC) -E -Xcompiler "-isystem $(CUDA_INSTALL_PATH)/include -MM" -o $@ $(CUDASRC)

include gpudirect.d

gpudirect.o: $(CUDASRC) gpudirect.d Makefile ../Makefile
	$(CUDACC) $(CUDAFLAGS) $(GLOBALVERSIONFLAGS) -c $< -o $@

clean:
	rm -f *.o

distclean:
	rm -f *.o

.PHONY: clean distclean gpudirect
