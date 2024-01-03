# A simple Makefile that passes all targets to subdirectories

# SUBDIRS := fast-cksum zeldovich-PLT
SUBDIRS := fast-cksum
MAKECMDGOALS ?= all

$(MAKECMDGOALS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(SUBDIRS)
