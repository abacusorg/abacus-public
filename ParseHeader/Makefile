CXX = icc
CXXFLAGS = -g
CPPFLAGS = -I../util
LEX = flex
YACC = bison

VPATH = ../util

# the lobrary name
LIBRARY = libparseheader.a

# these are the object files which will make up the object archive
CXXFILES =  phParser.tab.cc phScanner.cc ParseHeader.cc phDriver.cc HeaderStream.cc 
OFILES =  ${addprefix $(LIBRARY)(, ${addsuffix ), ${CXXFILES:.cc=.o}}}
HFILES = ParseHeader.hh  phDriver.hh  stringutil.hh phScanner.hh phParser.tab.hh

%.o:%.cc;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $*.o

(%.o):%.cc;
	$(CXX) ${CXFLAGS} $(CPPFLAGS) -c $< 
	$(AR) rv $@ $%
	- @ $(RM) -f $%

# the library
$(LIBRARY): $(OFILES) $(HFILES)
.PRECIOUS: $(LIBRARY)

test: 
	cd testdir && $(MAKE)

# the example executable

phexe: tryit.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^

HeaderStream.o: HeaderStream.cc ParseHeader.hh

ParseHeader.o: ParseHeader.cc ParseHeader.hh hDriver.hh

phParser.tab.o: phParser.tab.cc phParser.tab.hh

phScanner.o: phScanner.cc phParser.tab.hh phDriver.hh

phScanner.cc phScanner.hh: phScanner.ll phParser.tab.hh phDriver.hh
	$(LEX) $^

phParser.tab.cc phParser.tab.hh: phParser.yy phDriver.hh
	$(YACC) $<

phDriver.o: phDriver.cc phDriver.hh phScanner.hh phParser.tab.hh

phDriver.hh: phParser.tab.hh

ParseHeader.hh: phDriver.hh

.PHONY: clean distclean
clean:
	cd testdir && $(MAKE) clean
	$(RM) phDriver.o phParser.tab.o \
              phScanner.o phScanner.cc phScanner.hh \
              phParser.tab.cc phParser.tab.hh stack.hh location.hh position.hh \
              ParseHeader.o \
              tryit.o \
              phexe testit
	$(RM) $(wildcard *.o) $(wildcard *~) $(wildcard *.d)

distclean: 
	cd testdir && $(MAKE) clean
	$(RM) phDriver.o phParser.tab.o \
              phScanner.o phScanner.cc phScanner.hh \
              phParser.tab.cc phParser.tab.hh stack.hh location.hh position.hh \
              ParseHeader.o \
              tryit.o \
              phexe testit \
              libparseheader.a
	$(RM) $(wildcard *.o) $(wildcard *~) $(wildcard *.d)
