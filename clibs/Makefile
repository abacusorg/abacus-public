ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

LDFLAGS := $(CXXFLAGS) -shared $(PYTHON_LDFLAGS)
CXXFLAGS += -fPIC
CPPFLAGS += -I $(ROOT_DIR)/include/ -I $(ROOT_DIR)/Convolution $(THREAD_CPPFLAGS)

SO_TARGETS = libpermute.so liblightcones.so libreadabacus.so

all: libpermute.so libreadabacus.so

%.o:%.cpp
	@# Always use PYTHON_CXX for python extension compilations
	$(PYTHON_CXX) ${CXXFLAGS} $(CPPFLAGS) -MMD -c $<

$(SO_TARGETS): lib%.so: %.o
	$(PYTHON_CXX) $(CPPFLAGS) $(LDFLAGS) -o $@ $<

-include $(SO_TARGETS:.so=.d)

clean:
	$(RM) *.o *.so *.a
	
distclean: clean
	$(RM) *.d
	
.PHONY: clean distclean all
