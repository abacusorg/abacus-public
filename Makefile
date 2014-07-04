# Use 'True' or 'False' to toggle between GPU and CPU mode
USE_GPU='False'

#export CXX = icc -openmp -pthread -liomp5 -xHost -fp-model precise -fbuiltin -ip#-prof-use=weighted
export CXX = g++ -fopenmp -lgomp #-fprofile-use -fprofile-correction
GPUSPINFLAG = -DGPUTHREADFORCESPIN
#AVXFLAGS = -mavx -DAVXDIRECT -DAVXDIREC -DAVXMULTIPOLES
#GPUBLOCKINGFLAG = -DGPUBLOCKING
VERSIONFLAGS = -DDOUBLEPRECISION -DMAXCPD=8192 -DMAXSOURCELENGTH=1048576 $(GPUSPINFLAG) $(AVXFLAGS) $(GPUBLOCKINGFLAG)


ifeq ($(USE_GPU),'True')
# Use -DCUDADIRECT to use GPU; defaults to CPU otherwise
VERSIONFLAGS += -DCUDADIRECT
#VERSIONFLAGS += -DCUDA4
GPUDIRECT = gpudirect.o
endif
export VERSIONFLAGS

export CXXFLAGS= -O3 -DGITVERSION=\"`git rev-parse HEAD`\" $(VERSIONFLAGS) #-g #-DGLOBALPOS #-debug #-debug parallel 
# Could add -DGLOBALPOS here to switch the code to global positions.

CPPFLAGS = -I include -I Derivatives -I ParseHeader -I Library/include -I Library/lib/direct -I Library/lib/common -I/usr/local/cuda/include
CC_SRC = singlestep.cpp


-include ../Makefile.local
export ABACUS_VER = abacus_avx

LIBS =  -LParseHeader -LLibrary/lib -lparseheader -l$(ABACUS_VER) -lfftw3_omp -lfftw3 -ltbb
ifeq ($(USE_GPU),'True')
LIBS += gpudirect.o -L/usr/local/cuda/lib64  -lcudart -ltbb
endif

VPATH = singlestep : Convolution : Derivatives : python/clibs : zeldovich: Library/lib : Library/lib/direct

CLIBS = libpermute.so liblightcones.so

all: singlestep CreateDerivatives ConvolutionDriver zeldovich $(CLIBS) util tests powerspectrum

singlestep.o: singlestep.cpp lib$(ABACUS_VER).a Makefile $(GPUDIRECT)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -c -o $@ $<
	@sed -i 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\)\n\1 : \2,g' $*.d

gpudirect.o: gpu.cu Makefile
	nvcc --compiler-options $(GPUSPINFLAG) -DNFRADIUS=3 -I. -I/usr/local/cuda/include -ILibrary/include -arch compute_30 -code sm_30 -O3 -lineinfo -maxrregcount=48 -Xptxas="-v" -DUNIX -o $@ -c $<
	
libabacus_%.a:
	cd Library/lib && COMP=g++ _ABACUSDISTRIBUTION=$(ABACUS)/Library _ABACUSLIBRARY=libabacus_avx.a ./buildlibrary -static $(VERSIONFLAGS) -O3 -lfftw3


singlestep: singlestep.o $(GEN_OBJ) libparseheader.a lib$(ABACUS_VER).a Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o singlestep/$@ $< $(LIBS)

libparseheader.a:
	cd ParseHeader && $(MAKE) libparseheader.a
	
clean:
	#cd Library && $(MAKE) $@
	cd ParseHeader && $(MAKE) $@
	#cd Derivatives && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum &&$(MAKE) $@
	-$(RM) *.o *.d *.a *~

distclean:
	#cd Library && $(MAKE) $@
	cd ParseHeader && $(MAKE) $@
	#cd Derivatives && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum &&$(MAKE) $@
	-$(RM) *.o *.d *~ *.a a.out
	-$(RM) Library/lib/*.a


CreateDerivatives: lib$(ABACUS_VER).a
	cp Library/derivatives/* Derivatives/

ConvolutionDriver: convolutionwrapper.cpp include/Parameters.cpp
	cd Convolution && $(MAKE) $@
	
libpermute.so: perm.cpp
	cd python/clibs && $(MAKE) $@
	
liblightcones.so: lc.cpp
	cd python/clibs && $(MAKE) $@

util:
	cd util && $(MAKE) all

tests:
	cd Tests && $(MAKE) all	

powerspectrum:
	cd Analysis/PowerSpectrum &&$(MAKE)

zeldovich:zeldovich.cpp
	cd zeldovich && $(MAKE) $@
	
.PHONY: clean distclean generated_headers all zeldovich util tests powerspectrum 

-include $(CC_SRC:.cpp=.d)
