SUBDIRS = FoF PowerSpectrum Rockstar

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

%:
	for DIR in $(SUBDIRS); do \
		$(MAKE) -C $$DIR $@; \
	done
	