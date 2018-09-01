ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

all: Cartesian2Reduced.a $(MULTIPOLES_OBJ)

ETASM.o: generateCartesianTaylorASM.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -O0 -o genTaylorASM generateCartesianTaylorASM.c
	./genTaylorASM
	rm ./genTaylorASM
	$(CC) $(CPPFLAGS) $(CFLAGS) -c ETASM.c

CMASM.o: generateCartesianMultipolesASM.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -O0 -o genMultipoleASM generateCartesianMultipolesASM.c
	./genMultipoleASM
	rm genMultipoleASM
	$(CC) $(CPPFLAGS) $(CFLAGS) -c CMASM.c

CMAVX512.o: generateCartesianAVX512.py
	./generateCartesianAVX512.py
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c CMAVX512.cpp

ETAVX512.o: generateCartesianAVX512.py
	./generateCartesianAVX512.py
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c ETAVX512.cpp

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

EvaluateTaylorTest: EvaluateTaylor.cpp ETAVX512.o ETASM.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -xHost -Wall -fopenmp -DTEST $^ -o $@ -I../../include

EvaluateMultipolesTest: EvaluateMultipoles.cpp CMAVX512.o CMASM.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -xHost -Wall -fopenmp -DTEST $^ -o $@ -I../../include

clean:
	$(RM) ETASM.c CMASM.c CMAVX512.cpp *.o *.a
	$(RM) genTaylorASM genMultipoleASM CC2R
	$(RM) EvaluateTaylorTest EvaluateMultipolesTest
	$(RM) $(Cartesian2Reduced_OBJ) $(Cartesian2Reduced_SRC)

distclean: clean

.PHONY: all clean distclean

.INTERMEDIATE: $(Cartesian2Reduced_OBJ) $(Cartesian2Reduced_SRC)
