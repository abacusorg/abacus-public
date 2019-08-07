ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

all: Cartesian2Reduced.a $(MULTIPOLES_OBJ)

ETASM.o: generateCartesianTaylorASM.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -O0 -o genTaylorASM generateCartesianTaylorASM.c
	./genTaylorASM
	$(RM) ./genTaylorASM
	$(CC) $(CPPFLAGS) $(CFLAGS) -c ETASM.c

CMASM.o: generateCartesianMultipolesASM.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -O0 -o genMultipoleASM generateCartesianMultipolesASM.c
	./genMultipoleASM
	$(RM) genMultipoleASM
	$(CC) $(CPPFLAGS) $(CFLAGS) -c CMASM.c

CMAVX512%cpp ETAVX512%cpp: generateCartesianAVX512.py
	./generateCartesianAVX512.py

%AVX512.o: %AVX512.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

CM_unrolled%cpp ET_unrolled%cpp: generateCartesianUnrolled.py
	./generateCartesianUnrolled.py --maxorder=$(MAXORDER) --onlyorder=8

CM_VSX%cpp ET_VSX%cpp: generateCartesianVSX.py
	./generateCartesianVSX.py --maxorder=$(MAXORDER) --onlyorder=8

%_unrolled.o: %_unrolled.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fopt-info -Ofast -DUNROLLEDMULTIPOLES -Wall -c $<
%_VSX.o: %_VSX.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fopt-info -Ofast -Wall -c $<

Cartesian2Reduced_OBJ := $(shell echo Cartesian2Reduced{0..16}.o)
Cartesian2Reduced_TARGETPAT := $(Cartesian2Reduced_OBJ:.o=%cpp)
Cartesian2Reduced_SRC := $(Cartesian2Reduced_OBJ:.o=.cpp)

# Believe it or not, this tells make that all the targets are generated simultaneously!
# https://stackoverflow.com/a/3077254
$(Cartesian2Reduced_TARGETPAT): CreateCartesian2Reduced.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -O0 -o CC2R CreateCartesian2Reduced.cpp
	./CC2R
	$(RM) CC2R

Cartesian2Reduced.a: $(Cartesian2Reduced_OBJ)
	ar -r $@ $^

$(Cartesian2Reduced_OBJ): %.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(@:.o=.cpp)

EvaluateTaylorTest: EvaluateTaylor.cpp ETAVX512.o ETASM.o ET_unrolled.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Wall -fopenmp -DTEST $^ -o $@ -I../../include -lgsl -lgslcblas

EvaluateMultipolesTest: EvaluateMultipoles.cpp CMAVX512.o CMASM.o CM_unrolled.o CM_VSX.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Wall -fopenmp -DTEST $^ -o $@ -I../../include -lgsl -lgslcblas

clean:
	$(RM) ETASM.c CMASM.c ETAVX512.cpp CMAVX512.cpp CM_VSX.cpp CM_unrolled.cpp ET_unrolled.cpp *.o *.a
	$(RM) genTaylorASM genMultipoleASM CC2R
	$(RM) EvaluateTaylorTest EvaluateMultipolesTest 
	$(RM) $(Cartesian2Reduced_OBJ) $(Cartesian2Reduced_SRC)  CC2R

distclean: clean

.PHONY: all clean distclean

.INTERMEDIATE: $(Cartesian2Reduced_OBJ) $(Cartesian2Reduced_SRC) $(MULTIPOLES_OBJ:.o=.cpp)
