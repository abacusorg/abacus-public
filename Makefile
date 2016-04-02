-include Makefile.local

VPATH = singlestep : Convolution : Derivatives : python/clibs : zeldovich

CLIBS = libpermute.so liblightcones.so libreadabacus.so

all: singlestep CreateDerivatives ConvolutionDriver zeldovich $(CLIBS) util tests powerspectrum FoF

singlestep: libparseheader.a
	cd singlestep && $(MAKE) $@

CreateDerivatives:
	cd Derivatives && $(MAKE) $@

libparseheader.a:
	cd ParseHeader && $(MAKE) libparseheader.a
	
clean:
	cd ParseHeader && $(MAKE) $@
	cd Derivatives && $(MAKE) $@
	cd singlestep && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum && $(MAKE) $@
	cd Analysis/FoF && $(MAKE) $@
	-$(RM) *.o *.d *.a *~

distclean:	
	cd ParseHeader && $(MAKE) $@
	cd Derivatives && $(MAKE) $@
	cd singlestep && $(MAKE) $@
	cd Convolution && $(MAKE) $@
	cd python/clibs && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests/Spiral && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/PowerSpectrum &&$(MAKE) $@
	cd Analysis/FoF && $(MAKE) $@
	-$(RM) *.o *.d *~ *.a a.out



ConvolutionDriver: convolutionwrapper.cpp include/Parameters.cpp
	cd Convolution && $(MAKE) $@
	
libpermute.so: perm.cpp
	$(MAKE) -C python/clibs $@
	
liblightcones.so: lc.cpp
	$(MAKE) -C python/clibs $@
	
libreadabacus.so:
	$(MAKE) -C python/clibs $@

util:
	cd util && $(MAKE) all

tests:
	cd Tests && $(MAKE) all	
    
bench:
	cd singlestep && $(MAKE) $@

powerspectrum:
	cd Analysis/PowerSpectrum && $(MAKE)
    
FoF:
	cd Analysis/FoF && $(MAKE)

zeldovich: zeldovich.cpp
	cd zeldovich && $(MAKE) all
	
.PHONY: clean distclean generated_headers all zeldovich util tests powerspectrum singlestep

-include $(CC_SRC:.cpp=.d)
