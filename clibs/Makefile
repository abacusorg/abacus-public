CXX = g++
CXXFLAGS = -O3 -DDOUBLEPRECISION -fopenmp -shared -fPIC #-xhost
CPPFLAGS = -I../include/ -I../Convolution
LIBS = -lgomp #-liomp5


all: libpermute.so liblightcones.so libreadabacus.so

%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $<

libpermute.so: perm.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)
	
liblightcones.so: lc.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)
	
libreadabacus.so: read_abacus.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)
    
libreadabacus.a: read_abacus.o
	$(CXX) ${CXXFLAGS} -c -o $^ $(LIBS)
	ar rcs $@ $^

clean:
	rm -f *.o
	
distclean:
	rm -f *.o *.so *.a
	
.PHONY: clean distclean all
