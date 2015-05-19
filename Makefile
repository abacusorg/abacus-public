-include ../Makefile.local

VPATH = singlestep : Convolution : Derivatives : python/clibs : zeldovich

CLIBS = libpermute.so liblightcones.so

all: singlestep CreateDerivatives ConvolutionDriver zeldovich $(CLIBS) util tests powerspectrum

singlestep: libparseheader.a
	cd singlestep && $(MAKE) $@

CreateDerivatives:
	cd Derivatives && $(MAKE) $@

libparseheader.a:
	cd ParseHeader && $(MAKE) libparseheader.a
	
clean:
	cd ParseHeader && $(MAKE) $@
	cd singlestep && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum &&$(MAKE) $@
	-$(RM) *.o *.d *.a *~

distclean:	
	cd ParseHeader && $(MAKE) $@
	cd singlestep && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum &&$(MAKE) $@
	-$(RM) *.o *.d *~ *.a a.out



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
	
.PHONY: clean distclean generated_headers all zeldovich util tests powerspectrum singlestep

-include $(CC_SRC:.cpp=.d)
