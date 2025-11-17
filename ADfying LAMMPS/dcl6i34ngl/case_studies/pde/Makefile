SUBDIRS=a1s a1s_checkpointing cfd primal socfd t1s t2s_a1s t2s_a1s_checkpointing t2s_t1s

default: all

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all : $(SUBDIRS)
run : $(SUBDIRS)
check : $(SUBDIRS)
clean : $(SUBDIRS)
generate_results : $(SUBDIRS)

