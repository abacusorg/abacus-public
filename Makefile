VPATH = singlestep : Convolution : Derivatives : clibs : zeldovich

all: clibs singlestep CreateDerivatives ConvolutionDriver zeldovich util tests analysis AbacusCosmo

singlestep: libparseheader.a
	cd singlestep && $(MAKE) $@

CreateDerivatives:
	cd Derivatives && $(MAKE) $@

libparseheader.a:
	cd ParseHeader && $(MAKE) $@
	
clean:
	cd ParseHeader && $(MAKE) $@
	-cd Derivatives && $(MAKE) $@
	-cd singlestep && $(MAKE) $@
	-cd Convolution && $(MAKE) $@
	cd python/AbacusCosmo && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/ && $(MAKE) $@
	-$(RM) *.o *.d *.a *~

distclean:	
	cd ParseHeader && $(MAKE) $@
	-cd Derivatives && $(MAKE) $@
	-cd singlestep && $(MAKE) $@
	-cd Convolution && $(MAKE) $@
	cd python/AbacusCosmo && $(MAKE) $@
	cd zeldovich && $(MAKE) $@
	cd Tests && $(MAKE) $@
	cd util && $(MAKE) $@
	cd Analysis/ && $(MAKE) $@
	-$(RM) *.o *.d *~ *.a abacus.tar.gz
	-$(RM) singlestep/Makefile singlestep/Direct/Makefile singlestep/Multipoles/Makefile Convolution/Makefile Derivatives/Makefile Analysis/PowerSpectrum/Makefile Analysis/FoF/Makefile
	-$(RM) -rf autom4te.cache/ config.log config.status

ConvolutionDriver: convolutionwrapper.cpp include/Parameters.cpp
	cd Convolution && $(MAKE) $@
		
clibs:
	$(MAKE) -C clibs all

util:
	cd util && $(MAKE) all

tests:
	cd Tests && $(MAKE) all	
    
bench:
	cd singlestep && $(MAKE) $@

analysis:
	$(MAKE) -C Analysis/

zeldovich: zeldovich.cpp
	cd zeldovich && $(MAKE) all
    
AbacusCosmo:
	cd python/AbacusCosmo && $(MAKE) all

dist:
	$(RM) -rf .dist
	mkdir -p .dist/abacus
	cp -r * .dist/abacus
	$(MAKE) -C .dist/abacus distclean
	tar -C .dist -czf abacus.tar.gz --exclude='.*' abacus
	$(RM) -rf .dist
	
.PHONY:all clean distclean zeldovich util tests analysis singlestep dist AbacusCosmo clibs

-include $(CC_SRC:.cpp=.d)
