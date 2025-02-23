ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

all: Cartesian2Reduced.a $(MULTIPOLES_OBJ)

ETAVX.o: generateCartesianTaylorAVX.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -O0 -o genTaylorAVX generateCartesianTaylorAVX.cpp
	./genTaylorAVX
	$(RM) ./genTaylorAVX
	$(CC) $(CPPFLAGS) $(CFLAGS) -c ETAVX.c

CMAVX.o: generateCartesianMultipolesAVX.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -O0 -o genMultipoleAVX generateCartesianMultipolesAVX.cpp
	./genMultipoleAVX
	$(RM) genMultipoleAVX
	$(CC) $(CPPFLAGS) $(CFLAGS) -c CMAVX.c

CMAVX512%cpp ETAVX512%cpp: generateCartesianAVX512.py
	./generateCartesianAVX512.py --maxorder=$(MAXORDER)

%AVX512.o: %AVX512.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

CM_unrolled%cpp ET_unrolled%cpp: generateCartesianUnrolled.py
	./generateCartesianUnrolled.py --maxorder=$(MAXORDER)

CM_VSX%cpp ET_VSX%cpp: generateCartesianVSX.py
	./generateCartesianVSX.py --maxorder=$(MAXORDER)

%_unrolled.o: %_unrolled.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Ofast -Wall -DUNROLLEDMULTIPOLES -c $<
%_VSX.o: %_VSX.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Ofast -Wall -c $<

Cartesian2Reduced_OBJ := $(shell echo Cartesian2Reduced{0..$(MAXORDER)}.o Cartesian2ReducedDispatch.o)
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

EvaluateTaylorTest: EvaluateTaylor.cpp ETAVX512.o ETAVX.o ET_unrolled.o ET_VSX.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Wall -fopenmp -DTEST $^ -o $@ -I../../include $(GSL_LIBS)

EvaluateMultipolesTest: EvaluateMultipoles.cpp CMAVX512.o CMAVX.o CM_unrolled.o CM_VSX.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Wall -fopenmp -DTEST $^ -o $@ -I../../include $(GSL_LIBS)

clean:
	$(RM) *.o *.a
	$(RM) ETAVX.c CMAVX.c ETAVX512.cpp CMAVX512.cpp CM_VSX.cpp ET_VSX.cpp CM_unrolled.cpp ET_unrolled.cpp
	$(RM) genTaylorAVX genMultipoleAVX CC2R
	$(RM) EvaluateTaylorTest EvaluateMultipolesTest 
	$(RM) Cartesian2Reduced*.cpp

distclean: clean

.PHONY: all clean distclean

.INTERMEDIATE: $(Cartesian2Reduced_OBJ) $(Cartesian2Reduced_SRC) $(MULTIPOLES_OBJ:.o=.cpp)
