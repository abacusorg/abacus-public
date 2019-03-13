ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

CPPFLAGS += $(PARSEHEADER_CPPFLAGS) $(THREAD_CPPFLAGS) $(TCMALLOC_CPPFLAGS)
CPPFLAGS += -I $(ROOT_DIR)/Derivatives -I $(ROOT_DIR)/Convolution -I $(ROOT_DIR)/DirectIO\
	   		-I DataModel -I Direct -I IC -I Multipoles -I Output -I Timestep -I FOF

SINGLESTEP_LIBS += $(PARSEHEADER_LIBS) $(THREAD_LIBS)

MULTIPOLES_OBJ := $(addprefix Multipoles/, Cartesian2Reduced.a $(MULTIPOLES_OBJ))
EXTRA_OBJ = $(MULTIPOLES_OBJ)
    
ifeq ($(CUDA_ENABLED),CUDA_ENABLED)
	EXTRA_OBJ += Direct/gpudirect.a
	SINGLESTEP_LIBS += $(CUDA_LIBS)
	MAKE_DIRECT += _gpudirect
endif

PIPELINES = singlestep recover_multipoles benchmark_io standalone_fof
BOOTSTRAPS := $(addprefix _,$(PIPELINES))

all: _singlestep

-include $(addsuffix .d,$(PIPELINES))

# We need to run the recursion step (_gpudirect) first, before checking Direct/gpudirect.a
# But _gpudirect needs to be order-only w.r.t. singlestep, and unfortunately these seem to run last
# This fake parent enforces the desired ordering
# Warning: this means 'make singlestep' will not notice if the gpudirect.a has changed!
$(BOOTSTRAPS): _%: $(MAKE_DIRECT) %

$(PIPELINES): % : %.o $(EXTRA_OBJ)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(SINGLESTEP_LIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -c -o $@ $<

Direct/gpudirect.a: | _gpudirect
	$(MAKE) -C Direct $(notdir $@)

# This ensures we check if gpudirect.a needs to be remade without forcing a spurious re-link
_gpudirect:
	$(MAKE) -C Direct gpudirect.a

$(MULTIPOLES_OBJ): % :
	$(MAKE) -C Multipoles $(notdir $@)

clean: clean_recurse
	$(RM) *.o $(PIPELINES)

distclean: clean distclean_recurse
	$(RM) *.d

%_recurse:
	$(MAKE) -C Direct $*
	$(MAKE) -C Multipoles $*

.PHONY: all clean distclean cleanhere distcleanhere _gpudirect $(BOOTSTRAPS)
