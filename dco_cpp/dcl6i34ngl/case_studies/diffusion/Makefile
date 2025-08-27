SUBDIRS=obs f fd gt1s ga1s ga1s_sls sofd gt2s_gt1s gt2s_ga1s gt2s_ga1s_sls ga1s_cp ga1s_cp_sls gt2s_ga1s_cp gt2s_ga1s_cp_sls

default: all

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all : $(SUBDIRS)
run : $(SUBDIRS)
check : $(SUBDIRS)
clean : $(SUBDIRS)
generate_results : $(SUBDIRS)

