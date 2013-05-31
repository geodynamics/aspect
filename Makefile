# $Id: Makefile.large 22840 2010-11-22 22:28:16Z bangerth $

# The large projects Makefile looks much like the one for small
# projects. Basically, only the following seven parameters need to be
# set by you:

application-name  = aspect

# The next variable tells us the name of the executable. It is prefixed by
# `lib/' to designate its destination directory. Note that the program
# name depends on the dimension, so you can keep copies for the
# different dimensions around:
target   = lib/$(application-name)

# The `debug-mode' variable works as in the small projects Makefile:
debug-mode = on

# And so does the following variable. You will have to set it to
# something reasonable that, for example, includes the location where you
# put output files that you want the `make clean' rule to delete
clean-up-files =

# Finally, here is a variable which tells the `run' rule which
# parameters to pass to the executable. Usually, this will be the name
# of an input file.
run-parameters  = parameter-file.prm

# Now, this is the last variable you need to set, namely the path to
# the deal.II toplevel directory. DEAL_DIR is the old way of spelling
# this variable, please use DEAL_II_DIR instead. If nothing is specified 
# for either variable, we take ../../deal.II which works for some of
# us but is not likely to work for everyone who has a different directory
# layout
DEAL_DIR    ?= ../../deal.II
DEAL_II_DIR ?= $(DEAL_DIR)
D = $(DEAL_II_DIR)



#
#
# Usually, you will not need to change anything beyond this point.
#
#
# This tells `make' where to find the global settings and rules:
include $D/common/Make.global_options

CXXFLAGS.g += -Wunused-variable -Wextra

# list the directories and the various kinds of files
all-dirs := source \
	    source/simulator \
            source/termination_criteria \
            source/mesh_refinement \
            source/geometry_model \
            source/gravity_model \
            source/boundary_temperature \
            source/initial_conditions \
            source/compositional_initial_conditions \
            source/material_model \
            source/velocity_boundary_conditions \
            source/postprocess \
            source/postprocess/visualization

cc-files := $(shell for i in $(all-dirs) ; do echo $$i/*.cc ; done)

tmp1     := $(shell echo $(cc-files) | $(PERL) -pi -e 's,source/,,g; s,/,_,g;')
o-files  := $(addprefix lib/obj/, $(tmp1:.cc=.$(OBJEXT)) )
go-files := $(addprefix lib/obj/, $(tmp1:.cc=.g.$(OBJEXT)))

h-files     := $(wildcard include/aspect/*.h include/aspect/*/*h include/aspect/*/*/*h)
lib-h-files := $(shell echo $D/include/deal.II/*/*.h)

# As before, define two variables that denote the debug and optimized
# versions of the deal.II libraries:
libs.g   := $(lib-deal2.g)
libs.o   := $(lib-deal2.o)

INCLUDE += -Iinclude


# Now use the information from above to define the set of libraries to
# link with and the flags to be passed to the compiler:
ifeq ($(debug-mode),on)
  libraries = $(go-files) $(libs.g)
  flags     = $(CXXFLAGS.g)
else
  libraries = $(o-files) $(libs.o)
  flags     = $(CXXFLAGS.o)
endif


# The following two rules define how to compile C++ files into object
# files:
lib/obj/%.g.$(OBJEXT) :
	@echo =====$(application-name)=============debug=====$(MT)== $<
	@$(CXX) $(flags) -c $< -o $@
lib/obj/%.$(OBJEXT) :
	@echo =====$(application-name)=============optimized=$(MT)== $<
	@$(CXX) $(flags) -c $< -o $@



# Next define how to link the executable
build : $(target)$(EXEEXT)
$(target)$(EXEEXT) : $(libraries) Makefile
	@echo =====$(application-name)=======================$(MT)== Linking $@
	@$(CXX) -o $@ $(libraries) $(LIBS) $(LDFLAGS)



# Rule how to run the program
run: $(target)$(EXEEXT)
	./$(target)$(EXEEXT) $(run-parameters)

doc:
	@cd doc ; make

indent:
	@echo "============ Indenting all files"
	@astyle --options=lib/astyle.rc $(h-files) $(cc-files)

.PHONY: run build doc indent


# Rule how to clean up. This is split into several different rules to
# allow for parallel execution of commands:
clean: clean-lib clean-data
	-rm -f *~ */*~ */*/*~ source/Makefile.dep source/*/Makefile.dep
	-cd doc ; make clean
	-cd tests ; make clean

clean-lib:
	-rm -f lib/obj/*.$(OBJEXT) $(target)$(EXEEXT) lib/TAGS

clean-data:
	-rm -f $(clean-up-files)


# Again tell `make' which rules are not meant to produce files:
.PHONY: clean clean-data clean-lib run


# Rule to generate the dependency files, one for each source
# directory. These file are automagically remade whenever needed,
# i.e. whenever one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion below.
#
# If the command fails, then remove Makefile.dep again and fail
%/Makefile.dep: $(filter $(dir $@)%, $(cc-files)) \
                $(h-files) \
		$(lib-h-files) \
		Makefile
	@echo "====================================== Remaking $@"
	@(($D/common/scripts/make_dependencies -n $(INCLUDE) "-Blib/obj" \
			$(filter $(dir $@)%, $(cc-files)) \
		| $(PERL) -pe 's!debug/(.*)\.o:!$(subst source_,,$(subst /,_,$(dir $@)))\1.g.o:!g;' \
		| $(PERL) -pe 's!optimized/(.*)\.o:!$(subst source_,,$(subst /,_,$(dir $@)))\1.o:!g;' \
	  ) > $@) \
	 || (rm -f $@ ; false)


# include all the dependencies
include $(addsuffix /Makefile.dep, $(all-dirs))
