# A simple Makefile that passes all targets to subdirectories

# For external packages, we have to export variables since they
# don't include common.mk
ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

export CXX
export CXXFLAGS
export CPPFLAGS

SUBDIRS := fast-cksum
MAKECMDGOALS ?= all

$(MAKECMDGOALS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(SUBDIRS)
