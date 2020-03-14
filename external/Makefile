# A simple Makefile that passes all targets to subdirectories

SUBDIRS := trz fast-cksum zeldovich-PLT
MAKECMDGOALS ?= all

$(MAKECMDGOALS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(SUBDIRS)
