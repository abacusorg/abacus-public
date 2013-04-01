export CXX = icc
export CXXFLAGS = -mavx -O0 -g3 -openmp -DMAXCPD=8192 -DDOUBLE_PRECISION -DGITVERSION=\"`git rev-parse HEAD`\"

CPPFLAGS = -I Direct -I include -I Derivatives -I Multipoles -I Convolution -I ParseHeader -ILibrary/include
CC_SRC = singlestep.cpp


-include ../Makefile.local
ABACUS_VER = abacus_avx

LIBS = -LParseHeader -LLibrary -lparseheader -liomp5 -lfftw3 -l$(ABACUS_VER)

GEN_HDRS = externalmultipoles.h externaltaylor.h
GEN_OBJ = CMASM.o ETASM.o C2R.a

VPATH = singlestep : Direct : Multipoles : Convolution : Derivatives : python/clibs : zeldovich

CLIBS = libpermute.so liblightcones.so

all: singlestep CreateDerivatives ConvolutionDriver zeldovich $(CLIBS) util tests powerspectrum libabacus

singlestep: singlestep.o $(GEN_OBJ) libparseheader.a lib$(ABACUS_VER).a Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o singlestep/$@ $< $(addprefix Multipoles/,$(GEN_OBJ)) $(LIBS)


%.o: %.cpp | generated_headers Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -c -o $@ $<
	@sed -i 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\)\n\1 : \2,g' $*.d
	
CMASM.o: generateCartesianMultipolesASM.c
	cd Multipoles && $(MAKE) $@


ETASM.o: generateCartesianTaylorASM.c
	cd Multipoles && $(MAKE) $@


C2R.a: CreateCartesian2Reduced.cpp
	cd Multipoles && $(MAKE) $@
	
#$(GEN_HDRS):
#	cd Multipoles && $(MAKE) externaltaylor.h
	
generated_headers: $(GEN_HDRS)

libparseheader.a:
	cd ParseHeader && $(MAKE) libparseheader.a
	
clean:
	cd Multipoles && $(MAKE) $@
	cd ParseHeader && $(MAKE) $@
	cd Derivatives && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum &&$(MAKE) $@
	-$(RM) *.o *.d *~

distclean:
	cd Multipoles && $(MAKE) $@
	cd ParseHeader && $(MAKE) $@
	cd Derivatives && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum &&$(MAKE) $@
	-$(RM) *.o *.d *~ a.out


CreateDerivatives: CreateDerivatives.cpp
	cd Derivatives && $(MAKE) $@

ConvolutionDriver: ConvolutionDriver.cpp
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
