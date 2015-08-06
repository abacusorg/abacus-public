#export CXX = icc  -openmp-stubs -pthread -liomp5 #-xHost -fp-model precise -fbuiltin -ip#-prof-use=weighted
export CXX = g++ -fopenmp -lgomp --std=c++0x  #-fprofile-use -fprofile-correction 
export GLOBALVERSIONFLAGS = -DFLOATPRECISION -DNFRADIUS=2 -DDIRECTSPLINE
export VERSIONFLAGS = $(GLOBALVERSIONFLAGS) -mavx -DAVXMULTIPOLES\
	       -DMAXCPD=8192 -DMAXSOURCELENGTH=1048576 -O3 -DIOTHREADED -DCUDADIRECT #-DAVXDIRECT #-ggdb3
export CXXFLAGS= -DGITVERSION=\"`git rev-parse HEAD`\" $(VERSIONFLAGS) -Wformat=0
# Could add -DGLOBALPOS here to switch the code to global positions.
# Use -DDIRECTSPLINE for spline softening

-include ../Makefile.local

CPPFLAGS = -I ../include -I ../Derivatives -I ../ParseHeader -I../Convolution\
	   -I DataModel -I Direct/ -I IC -I Multipoles -I Output -I Timestep\
	   -I/usr/local/cuda/include


LIBS =  -L../ParseHeader -lparseheader -lfftw3_omp -lfftw3 -ltbb\
	Multipoles/ETASM.o Multipoles/CMASM.o Multipoles/Cartesian2Reduced.a
    
ifneq ($(DISABLE_GPU),1)
LIBS += Direct/gpudirect.o -L/usr/local/cuda/lib64 -lcudart
endif


GEN_OBJ = ETASM.o CMASM.o Cartesian2Reduced.a

VPATH = Direct : Multipoles

singlestep: singlestep.o $(GEN_OBJ) gpudirect.a | singlestep.d 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS)

singlestep.d: singlestep.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -M -MF $@ $<

include singlestep.d

singlestep.o: singlestep.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c -o $@ $<

gpudirect.a:
	cd Direct && $(MAKE)

$(GEN_OBJ):
	cd Multipoles && $(MAKE) $@


clean:
	cd Direct && $(MAKE) $@
	-$(RM) *.o *.d

distclean:
	cd Direct && $(MAKE) $@
	-$(RM) *.o *.d singlestep


.PHONY: clean distclean gpudirect.o
