all: spiral

RECURSIVE = clean distclean

spiral:
	$(MAKE) -C Spiral
	
$(RECURSIVE): %:
	$(MAKE) -C Spiral $@

.PHONY: spiral clean distclean all
