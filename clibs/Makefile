ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

LDFLAGS := $(CXXFLAGS) -shared $(PYTHON_LDFLAGS)
CXXFLAGS += -fPIC
PYTHON_CPPFLAGS += -I $(ROOT_DIR)/DirectIO $(THREAD_CPPFLAGS) -O3

SO_TARGETS = libpermute.so liblightcones.so libreadabacus.so

all: libpermute.so libreadabacus.so

# Always use PYTHON_CXX for python extension compilations
%.o:%.cpp
	$(PYTHON_CXX) ${CXXFLAGS} $(PYTHON_CPPFLAGS) -MMD -c $<

$(SO_TARGETS): lib%.so: %.o
	$(PYTHON_CXX) $(PYTHON_CPPFLAGS) $(LDFLAGS) -o $@ $<

-include $(SO_TARGETS:.so=.d)

clean:
	$(RM) *.o *.so *.a
	
distclean: clean
	$(RM) *.d
	
.PHONY: clean distclean all
