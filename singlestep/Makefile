ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

CPPFLAGS += $(PARSEHEADER_CPPFLAGS) $(THREAD_CPPFLAGS) -I $(ROOT_DIR)/Derivatives -I $(ROOT_DIR)/Convolution\
	   -I DataModel -I Direct -I IC -I Multipoles -I Output -I Timestep -I FOF

LIBS += $(PARSEHEADER_LIBS) $(THREAD_LIBS)

MULTIPOLES_OBJ := $(addprefix Multipoles/, Cartesian2Reduced.a $(MULTIPOLES_OBJ))
EXTRA_OBJ = $(MULTIPOLES_OBJ)
    
ifeq ($(CUDA_ENABLED),CUDA_ENABLED)
	EXTRA_OBJ += Direct/gpudirect.a
	LIBS += $(CUDA_LIBS)
endif

PIPELINES = singlestep make_multipoles benchmark_io standalone_fof

all: singlestep

-include $(addsuffix .d,$(PIPELINES))

$(PIPELINES): % : %.o $(EXTRA_OBJ)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -MMD $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c -o $@ $<

Direct/gpudirect.a:
	$(MAKE) -C Direct $(notdir $@)

$(MULTIPOLES_OBJ): % :
	$(MAKE) -C Multipoles $(notdir $@)

clean: clean_recurse
	$(RM) *.o $(PIPELINES)

distclean: clean distclean_recurse
	$(RM) *.d

%_recurse:
	$(MAKE) -C Direct $*
	$(MAKE) -C Multipoles $*

.PHONY: all clean distclean cleanhere distcleanhere
