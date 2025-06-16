SUBDIRS = PowerSpectrum

.PHONY: $(SUBDIRS) all default

# don't build Rockstar by default (expensive)
default: PowerSpectrum

$(SUBDIRS):
	$(MAKE) -C $@
	
all: $(SUBDIRS)

# Pass on all targets to children
%:
	for DIR in $(SUBDIRS); do \
		$(MAKE) -C $$DIR $@; \
	done
