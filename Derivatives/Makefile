CXX= g++ -O3 -fopenmp -I../include -I../singlestep/Multipoles/ 

CreateDerivatives:CreateDerivatives.cpp
	$(CXX) $< -o $@ -lfftw3 -lm -lgomp

clean:
	$(RM) ./CreateDerivatives
    
distclean: clean