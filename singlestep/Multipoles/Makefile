CXX=g++ -O1 -DAVXMULTIPOLES
INCLUDES=-I../../include

ETASM.o: cartesian2reduced.cpp
	$(CXX) $(INCLUDES) -o genTaylorASM generateCartesianTaylorASM.c
	./genTaylorASM
	rm ./genTaylorASM
	$(CXX) $(INCLUDES) -c ETASM.c

CMASM.o: cartesian2reduced.cpp
	$(CXX) $(INCLUDES) -o genMultipoleASM generateCartesianMultipolesASM.c
	./genMultipoleASM
	rm genMultipoleASM
	$(CXX) $(INCLUDES) -c CMASM.c

Cartesian2Reduced.a: CreateCartesian2Reduced.cpp
	$(CXX) $(INCLUDES) -o CC2R CreateCartesian2Reduced.cpp
	./CC2R
	rm CC2R
	ar -r $@ *2Reduced*.o

clean:
	rm -f ETASM.c CMASM.c *.o 

distclean: clean

.PHONY:clean distclean
