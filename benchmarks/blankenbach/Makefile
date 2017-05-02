# Makefile to run all Blankenbach benchmarks
#
# This Makefile first generates caseXX_refY.prm files based on the existing
# base_caseXX.prm files for the specified global refinement levels Y (see
# "refinements" variable below).
#
# In a second step, all simulations for caseXX_refY.prm files are run and the
# statistics (that include RMS and Nusselt numbers) are saved as
# caseXX_refY.stat.
#
# Finally, using the helper script "get_statistics", the relevant columns are
# extracted and printed to the screen.
#
# Note: you can use something like "make -j 4" to run more than one
# computation at once or you can run a specific computation via "make
# case1a_ref5.stat".


# specify all refinements to run here:
refinements:=4 5 6 7

# bases_without is "case1a case1b ..."
bases:=$(wildcard base_case*.prm)
bases_without:=$(subst base_,,$(basename $(bases)))

# create a list of "caseXX_refY.stat" names:
files:=$(subst base_,,$(bases))
prmfiles:=$(foreach ref,$(refinements),$(subst .prm,_ref$(ref).prm,$(files)))
statfiles:=$(addsuffix .stat, $(basename $(prmfiles)))

# number of processors to run with
nproc=1

# main target to show the statistics
all: $(statfiles)
	./get_statistics *.stat

# keep the prm files around (otherwise make will delete intermediate files)
.PRECIOUS: $(foreach ref,$(refinements),%_ref$(ref).prm)

# The rule to generate caseXX_refY.prm
$(foreach ref,$(refinements),%_ref$(ref).prm): base_%.prm
	@for ref in $(refinements) ; do \
		f="$*_ref$$ref.prm";\
		echo "making $$f";\
		echo "include base_$*.prm" > $$f; \
		echo "set Additional shared libraries= ./plugin/libblankenbach.so" >>$$f; \
		echo "  subsection Mesh refinement" >> $$f; \
		echo "set Initial global refinement=$$ref" >> $$f; \
		echo "end" >> $$f; \
		echo "set Output directory=output-$*_ref$${ref}" >> $$f; \
	done 

# The rule to run ASPECT to generate caseXX_refY.stat
%.stat: %.prm
	@echo "** $*.prm -> output-$*/statistics -> $@"
	mpirun -n $(nproc) ./plugin/aspect $*.prm >/dev/null
	@cp output-$*/statistics $@
	#@rm -rf output-$*
	@echo created $@.

.PHONY: all
