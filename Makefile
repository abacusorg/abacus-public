export CXX = icc -openmp -liomp5
#export CXX = gcc -fopenmp -lgomp
export VERSIONFLAGS = -DFLOATPRECISION -DAVXDIRECT -DAVXMULTIPOLES

export CXXFLAGS = -mavx -O0 -g3 -DMAXCPD=8192 -mavx -DMAXCPD=8192 -DMAXSOURCELENGTH=1048576 -DGITVERSION=\"`git rev-parse HEAD`\" $(VERSIONFLAGS)
# Could add -DGLOBALPOS here to switch the code to global positions.

CPPFLAGS = -I include -I Derivatives -I ParseHeader -ILibrary/include
CC_SRC = singlestep.cpp


-include ../Makefile.local
export ABACUS_VER = abacus_avx

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
	-$(RM) *.o *.d *.a *~

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
	-$(RM) *.o *.d *~ *.a a.out


CreateDerivatives: CreateDerivatives.cpp
	cd Derivatives && $(MAKE) $@

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
