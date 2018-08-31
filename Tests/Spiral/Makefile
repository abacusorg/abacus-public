ROOT_DIR := ../..
include $(ROOT_DIR)/common.mk

CXXFLAGS := -O1 $(THREAD_LIBS)
CPPFLAGS += $(THREAD_CPPFLAGS) -DDOUBLE_PRECISION
LIBS := -lfftw3  # in theory this could be MKL, but we won't think that hard

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
