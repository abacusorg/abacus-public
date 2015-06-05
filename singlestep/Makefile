#export CXX = icc  -openmp-stubs -pthread -liomp5 #-xHost -fp-model precise -fbuiltin -ip#-prof-use=weighted
export CXX = g++ -fopenmp -lgomp --std=c++0x  #-fprofile-use -fprofile-correction 
export GLOBALVERSIONFLAGS = -DFLOATPRECISION -DNFRADIUS=3 #-DDIRECTSPLINE
export VERSIONFLAGS = $(GLOBALVERSIONFLAGS) -mavx -DAVXMULTIPOLES\
	   -DCUDADIRECT -DAVXDIRECT -DMAXCPD=8192 -DMAXSOURCELENGTH=1048576 -O0 -ggdb3
export CXXFLAGS= -DGITVERSION=\"`git rev-parse HEAD`\" $(VERSIONFLAGS) -Wformat=0
# Could add -DGLOBALPOS here to switch the code to global positions.
# Use -DDIRECTSPLINE for spline softening


CPPFLAGS = -I ../include -I ../Derivatives -I ../ParseHeader -I../Convolution\
	   -I DataModel -I Direct/ -I IC -I Multipoles -I Output -I Timestep\
	   -I/usr/local/cuda/include


LIBS =  -L../ParseHeader -lparseheader -lfftw3_omp -lfftw3 -ltbb\
	Direct/gpudirect.o Multipoles/ETASM.o Multipoles/CMASM.o Multipoles/Cartesian2Reduced.a\
	-L/usr/local/cuda/lib64 -lcudart -ltbb

GEN_OBJ = ETASM.o CMASM.o Cartesian2Reduced.a
-include ../Makefile.local

VPATH = Direct : Multipoles

singlestep: singlestep.o $(GEN_OBJ) gpudirect.o | singlestep.d 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS)
	$(RM) singlestep.d

singlestep.d:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -M -MF $@ $<

-include singlestep.d

singlestep.o: singlestep.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c -o $@ $<

gpudirect.o:
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
