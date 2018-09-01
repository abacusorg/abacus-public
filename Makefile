-include common.mk

all: clibs singlestep CreateDerivatives ConvolutionDriver zeldovich util tests analysis AbacusCosmo

singlestep: ParseHeader
	$(MAKE) -C singlestep all

CreateDerivatives:
	$(MAKE) -C Derivatives $@
	
clean: clean_recurse

distclean: clean distclean_recurse
	$(RM) abacus.tar.gz
	$(RM) common.mk
	$(RM) -rf autom4te.cache/ config.log config.status

# One annoyance with common.mk is that we don't want Makefiles to build without it,
# but we would prefer they be able to run clean.
# We could require a second, shadow common.mk be included that prevents CXX from running,
# but that's another line of text in each Makefile and isn't very clean
# So instead we'll just accept that they'll fail

%_recurse:
	-$(MAKE) -C ParseHeader $*
	-$(MAKE) -C Derivatives $*
	-$(MAKE) -C singlestep $*
	-$(MAKE) -C Convolution $*
	-$(MAKE) -C Abacus/Cosmology $*
	-$(MAKE) -C zeldovich-PLT $*
	-$(MAKE) -C Tests $*
	-$(MAKE) -C util $*
	-$(MAKE) -C Analysis $*
	-$(MAKE) -C clibs $*

ConvolutionDriver: ParseHeader
	$(MAKE) -C Convolution $@
		
clibs:
	$(MAKE) -C clibs all

util:
	$(MAKE) -C util all

tests:
	$(MAKE) -C Tests all	

analysis: clibs ParseHeader util AbacusCosmo
	$(MAKE) -C Analysis

zeldovich:
	$(MAKE) -C zeldovich-PLT all

AbacusCosmo:
	$(MAKE) -C Abacus/Cosmology all

# Make an abacus.tar.gz file for distribution/archival purposes
dist:
	$(RM) -rf .dist
	mkdir -p .dist/abacus
	cp -r * .dist/abacus
	$(MAKE) -C .dist/abacus distclean
	tar -C .dist -czf abacus.tar.gz --exclude='.*' abacus
	$(RM) -rf .dist

# Usually you do not need to invoke the following directly; it is used by ./configure
tcmalloc: gperftools/lib/libtcmalloc_minimal.so
gperftools/lib/libtcmalloc_minimal.so:
	@echo "Building tcmalloc... this may take a minute but will only be done once"
	@cd gperftools && \
	./configure --enable-minimal --prefix=$(shell pwd)/gperftools > /dev/null && \
	make > /dev/null && make install > /dev/null
	
.PHONY:all clean distclean zeldovich util tests analysis singlestep dist AbacusCosmo clibs ConvolutionDriver tcmalloc ParseHeader
