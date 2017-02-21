# Default: just make PowerSpectrum
SUBDIRS = PowerSpectrum FoF Rockstar

.PHONY: $(SUBDIRS) all

$(SUBDIRS):
	$(MAKE) -C $@
	
all: $(SUBDIRS)

%:
	for DIR in $(SUBDIRS); do \
		$(MAKE) -C $$DIR $@; \
	done
	