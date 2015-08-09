ifneq ($(DISABLE_GPU),1)

CUDACC=nvcc
#CUDAFLAGS= -I/usr/local/cuda/include -m64\
	    --generate-code arch=compute_35,code=sm_35\
	    --generate-code arch=compute_30,code=sm_30\
	    -lineinfo
CUDASRC= gpudirect.cu

Blocksize = 64
Tile = 3
ILP = 1
Unroll = 1
Tilesize = $(shell echo $(Blocksize)\*$(Tile) | bc)

OPTIONS=-DDirectNT=$(Blocksize) -DDirectCPT=$(Tile) #-DPTXDIRECT
KERNELPTXFLAGS=-v -O1 --compile-only  -maxrregcount=54 --generate-line-info
PTXFLAGS= -v -O1 --compile-only  -maxrregcount=54 --generate-line-info
CUDASM= sm_35
CUDAFLAGS= $(OPTIONS) $(GLOBALVERSIONFLAGS) -m64 -arch compute_35 -code $(CUDASM) -Xptxas="$(PTXFLAGS)" -Xcompiler -fpic,-O3,-march=native,-mavx

gpudirect: gpudirect.a gpudirect.d bench

gpudirect.d:
	$(CUDACC) -E -Xcompiler "-isystem $(CUDA_INSTALL_PATH)/include -MM" -o $@ $(CUDASRC)

include gpudirect.d

GeneratedASMKernel.ptx: ./GenPTXKernel.py Makefile
	./GenPTXKernel.py $(Tilesize) $(ILP) $(Unroll)

GeneratedASMKernel.cubin: GeneratedASMKernel.ptx
	ptxas $(KERNELPTXFLAGS) --gpu-name $(CUDASM) $< -o $@

ASMDirectTile.o: GeneratedASMKernel.cubin
	export GAK_FH="`filehash -s "$(PTXFLAGS)" "$(<:.o=.ptx)"`"
	mkdir -p nvtemp
	nvcc --keep --keep-dir nvtemp --device-c $(CUDAFLAGS) -DDUMMYPTX $(@:.o=.cu) -o $@
	fatbinary --create="$(@:.o=.fatbin)" -64 \
                --key="$(GAK_FH)" --cmdline="$(PTXFLAGS)"\
                "--image=profile=$(CUDASM),file=$<" --cuda --device-c\
                --embedded-fatbin="$(@:.o=.fatbin.c)"
	mv $(@:.o=.fatbin.c) nvtemp/ASMDirectTile.fatbin.c
	gcc -D__CUDA_ARCH__=350 -E -x c++\
                -DCUDA_DOUBLE_MATH_FUNCTIONS -D__CUDA_PREC_DIV -D__CUDA_PREC_SQRT\
                "-I/usr/local/cuda/bin/..//include"   -m64 "nvtemp/ASMDirectTile.cudafe1.cpp"\
                > "nvtemp/ASMDirectTile.cu.cpp.ii"
	gcc -c -x c++ "-I/usr/local/cuda/bin/..//include"\
                -fpreprocessed -m64 -o $@ "nvtemp/ASMDirectTile.cu.cpp.ii"


gpudirect.dev.o: $(CUDASRC) gpudirect.d Makefile ../Makefile
	$(CUDACC) $(CUDAFLAGS) --device-c $< -o $@

devkernel.o: gpudirect.dev.o ASMDirectTile.o
	$(CUDACC) $(CUDAFLAGS) -dlink $^ -o $@

gpudirect.a: devkernel.o gpudirect.dev.o ASMDirectTile.o
	$(CUDACC) $(CUDAFLAGS) -lib $^ -o $@

bench: DirectBench.cpp gpudirect.a
	$(CUDACC) $(CUDAFLAGS) -I$(ABACUS)/include $^ -o DirectBench -lm -lgsl -lgslcblas -Xcompiler -fopenmp
	#./DirectBench 147 42 5 1
	#./DirectBench 275 125 40 2

endif

clean:
	rm -f *.o DirectBench

distclean:
	rm -f *.o *.d *.a DirectBench

.PHONY: clean distclean gpudirect
