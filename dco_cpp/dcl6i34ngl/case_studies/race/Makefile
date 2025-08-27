SUBDIRS=cfd ga1s ga1v ga2s_ga1s ga2s_gt1s gt1s gt1v gt2s_ga1s gt2s_gt1s socfd

default: all

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all : $(SUBDIRS)
run : $(SUBDIRS)
check : $(SUBDIRS)
clean : $(SUBDIRS)
generate_results : $(SUBDIRS)

