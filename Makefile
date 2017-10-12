VPATH = singlestep : Convolution : Derivatives : clibs : zeldovich-PLT

all: clibs singlestep CreateDerivatives ConvolutionDriver zeldovich util tests analysis AbacusCosmo

singlestep: libparseheader.a
	cd singlestep && $(MAKE) all

CreateDerivatives:
	cd Derivatives && $(MAKE) $@

libparseheader.a:
	cd ParseHeader && $(MAKE) $@
	
clean:
	cd ParseHeader && $(MAKE) $@
	-cd Derivatives && $(MAKE) $@
	-cd singlestep && $(MAKE) $@
	-cd Convolution && $(MAKE) $@
	cd python/Abacus/Cosmology && $(MAKE) $@
	cd zeldovich-PLT && $(MAKE) $@
	cd Tests && $(MAKE) $@
	cd util && $(MAKE) $@
	-cd Analysis/ && $(MAKE) $@
	-cd clibs/ && $(MAKE) $@
	-$(RM) *.o *.d *.a *~

distclean:	
	cd ParseHeader && $(MAKE) $@
	-cd Derivatives && $(MAKE) $@
	-cd singlestep && $(MAKE) $@
	-cd Convolution && $(MAKE) $@
	cd python/Abacus/Cosmology && $(MAKE) $@
	cd zeldovich-PLT && $(MAKE) $@
	cd Tests && $(MAKE) $@
	cd util && $(MAKE) $@
	-cd Analysis/ && $(MAKE) $@
	-cd clibs/ && $(MAKE) $@
	-$(RM) *.o *.d *~ *.a abacus.tar.gz
	-$(RM) singlestep/Makefile singlestep/Direct/Makefile singlestep/Multipoles/Makefile Convolution/Makefile Derivatives/Makefile Analysis/PowerSpectrum/Makefile Analysis/FoF/Makefile clibs/Makefile
	-$(RM) -rf autom4te.cache/ config.log config.status

ConvolutionDriver: libparseheader.a
	cd Convolution && $(MAKE) $@
		
clibs:
	$(MAKE) -C clibs all

util:
	cd util && $(MAKE) all

tests:
	cd Tests && $(MAKE) all	
    
bench:
	cd singlestep && $(MAKE) $@

analysis: clibs libparseheader.a util AbacusCosmo
	$(MAKE) -C Analysis/

zeldovich: zeldovich.cpp
	cd zeldovich-PLT && $(MAKE) all
    
AbacusCosmo:
	cd python/Abacus/Cosmology && $(MAKE) all

dist:
	$(RM) -rf .dist
	mkdir -p .dist/abacus
	cp -r * .dist/abacus
	$(MAKE) -C .dist/abacus distclean
	tar -C .dist -czf abacus.tar.gz --exclude='.*' abacus
	$(RM) -rf .dist
	
.PHONY:all clean distclean zeldovich util tests analysis singlestep dist AbacusCosmo clibs ConvolutionDriver

-include $(CC_SRC:.cpp=.d)
