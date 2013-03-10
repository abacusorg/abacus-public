#CXX = icc
CXXFLAGS = -O1
INCLUDES = -I../../include/
LIBS = 

CFLAGS = -O1 ${INCLUDES}

TARGET = phheader 

all: ${TARGET}

%.o:%.cpp;
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $< $(INCLUDES)

clean:
	rm -f *.o
	
distclean:
	rm -f *.o ${TARGET}
	
.PHONY: clean distclean all
