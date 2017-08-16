SUBDIRS = PowerSpectrum FoF 

.PHONY: $(SUBDIRS) all default

# don't build Rockstar by default (expensive)
default: PowerSpectrum FoF

$(SUBDIRS):
	$(MAKE) -C $@
	
all: $(SUBDIRS)

%:
	for DIR in $(SUBDIRS); do \
		$(MAKE) -C $$DIR $@; \
	done
	
