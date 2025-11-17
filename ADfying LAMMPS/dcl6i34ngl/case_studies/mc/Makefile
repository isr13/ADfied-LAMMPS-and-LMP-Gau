SUBDIRS=a1s a1s_ensemble a1s_preaccumulator cfd primal socfd t1s t2s_a1s t2s_a1s_ensemble t2s_t1s

default: all

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all : $(SUBDIRS)
run : $(SUBDIRS)
check : $(SUBDIRS)
clean : $(SUBDIRS)
generate_results : $(SUBDIRS)

