ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

CXXFLAGS += $(THREAD_LIBS)
CPPFLAGS += $(THREAD_CPPFLAGS) -DDOUBLE_PRECISION
LIBS := $(FFTW_LIBS) -lstdc++fs

all:makespiralics makeanalytic

makespiralics: spiralic.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)

makeanalytic: AnalyticSpiral.o
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -c $<
	
clean:
	$(RM) *.o makespiralics makeanalytic analytic

distclean: clean

.PHONY: clean distclean all
