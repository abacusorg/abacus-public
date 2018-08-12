ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

LDFLAGS := $(CXXFLAGS) -shared -Wl,--as-needed
CXXFLAGS += -fPIC
CPPFLAGS += -I../include/ -I../Convolution

SO_TARGETS = libpermute.so liblightcones.so libreadabacus.so

all: libpermute.so libreadabacus.so

%.o:%.cpp
	$(CXX) ${CXXFLAGS} $(CPPFLAGS) -MMD -c $<

$(SO_TARGETS): lib%.so: %.o
	$(CXX) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

-include $(SO_TARGETS:.so=.d)

clean:
	$(RM) *.o *.so *.a
	
distclean: clean
	$(RM) *.d
	
.PHONY: clean distclean all
