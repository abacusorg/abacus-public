export CXX = icc
export CXXFLAGS = -mavx -O0 -g3 -openmp -DMAXCPD=8192 -DDOUBLE_PRECISION -mavx -DAVXDIRECT -DAVXMULTIPOLES -DMAXCPD=8192 -DMAXSOURCELENGTH=1048576 -DGITVERSION=\"`git rev-parse HEAD`\"

CPPFLAGS = -I Direct -I include -I Derivatives -I Multipoles -I ParseHeader -ILibrary/include
CC_SRC = singlestep.cpp


-include ../Makefile.local
ABACUS_VER = abacus_avx

LIBS = -LParseHeader -LLibrary -lparseheader -liomp5 -lfftw3 $(ABACUS_VER).a


VPATH = singlestep : Convolution : Derivatives : python/clibs : zeldovich

CLIBS = libpermute.so liblightcones.so

all: singlestep CreateDerivatives ConvolutionDriver zeldovich $(CLIBS) util tests powerspectrum

singlestep: singlestep.o $(GEN_OBJ) libparseheader.a $(ABACUS_VER).a Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o singlestep/$@ $< $(LIBS)


%.o: %.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -c -o $@ $<
	@sed -i 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\)\n\1 : \2,g' $*.d

abacus_%.a:
	cd Library && $(MAKE) $@
	

libparseheader.a:
	cd ParseHeader && $(MAKE) libparseheader.a
	
clean:

	cd Library && $(MAKE) $@
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
	cd Library && $(MAKE) $@
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
