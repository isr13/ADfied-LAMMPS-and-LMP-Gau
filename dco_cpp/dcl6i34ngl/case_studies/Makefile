SUBDIRS=diffusion pde mc race

default: all

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all : $(SUBDIRS)
run : $(SUBDIRS)
check : $(SUBDIRS)
clean : $(SUBDIRS)
generate_results : $(SUBDIRS)

