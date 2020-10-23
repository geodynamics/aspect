#
# This Makefile compiles all plugins and runs all prm files in the
# subdirectories of the benchmarks folder
#
# usage:
#   make -f check.mk BUILD=/ssd/aspect-git/build -j4
#
# optionally specify "CHECK=--validate" to only check the .prm's and not run them.

SHELL := /bin/bash

def= \
if [[ ! -e ${BUILD}/AspectConfig.cmake || ! -e ${BUILD}/aspect  ]]; \
then \
    echo "please run using     make -f check.mk BUILD=/aspect/build/directory -j4"; \
    exit 1; \
fi; \
make_lib() { \
cd $$1; \
if [[ -e CMakeLists.txt ]]; \
then \
  echo "building plugin in `pwd` using ${BUILD}..."; \
  rm -rf CMakeCache.txt CMakeFiles; \
  cmake -D Aspect_DIR=${BUILD} -G "Unix Makefiles" -D CMAKE_CXX_FLAGS='-Werror' . >/dev/null || { echo "cmake in `pwd` failed!"; return 1; }; \
  make >/dev/null || { echo "make in `pwd` failed!"; return 2; }; \
  echo "done building plugin in `pwd`"; \
fi; \
};\
run_prm() { \
cd $$1; \
prm=$$2; \
echo "    running $$prm in `pwd` ..."; \
[[ -f $$prm ]] || { echo "File '$$prm' missing in `pwd`."; return 7; }; \
cp $$prm $$prm.tmp || return 5; \
echo "set End time=0" >> $$prm.tmp; \
echo "set Max nonlinear iterations = 5" >> $$prm.tmp; \
${BUILD}/aspect ${CHECK} $$prm.tmp >/dev/null || { rm -f $$prm.tmp; echo "run of $$prm failed"; return 2; }; \
rm -f $$prm.tmp; \
}; \
run_all_prms() { \
cd $$1; \
echo "  running all prms in `pwd` ..."; \
for prm in *.prm; \
  do \
    run_prm . $$prm || return 4; \
  done; \
echo "  running all prms in `pwd` done"; \
}

allsubdirs:= $(wildcard */)
subdirs:= $(filter-out output-%,$(allsubdirs))

all: main $(subdirs)

# You can not mix implicit rules (like %/ below) with phony targets, instead,
# define a phony "dummy" target that all other targets depend on. This will
# always rebuild them.
.PHONY: all dummy
dummy:

# default rule for each folder: compile lib and run all prms
%/: dummy
	+@$(def); make_lib $@
	@$(def); run_all_prms $@

# there are no .prm files in the main directory:
main: dummy
	@echo "ok"



# custom rules. Make them dependent on dummy like this:
#
# example/: dummy
#	@$(def); run_prm $@ test.prm

blankenbach/: dummy
	+@$(def); make_lib $@/plugin
	@$(def); run_all_prms $@

crameri_et_al/:  dummy
	+@$(def); make_lib $@/case_1
	@$(def); run_all_prms $@/case_1
	@$(def); run_all_prms $@/case_2

davies_et_al/: dummy
	+@$(def); make_lib $@/case-2.3-plugin
	@$(def); run_prm $@ case-2.1.prm
	@$(def); run_all_prms $@

entropy_adiabat/: dummy
	+@$(def); make_lib $@/plugins
	@$(def); run_all_prms $@

free_surface_tractions/: dummy
	@$(def); run_all_prms $@/viscoelastic
	@$(def); run_all_prms $@/viscous

# default file is too fine:
geoid-spectral-comparison/: dummy
	@$(def); run_prm $@ spectral-comparison.prm "subsection Mesh refinement\n set Initial global refinement=2\n end"

inclusion/: dummy
	+@$(def); make_lib $@
	@$(def); run_prm $@ global.prm.base
	@$(def); run_prm $@ adaptive.prm.base
	+@$(def); make_lib $@/compositional_fields
	@$(def); run_all_prms $@/compositional_fields


# TODO: no .prm in folder:
nonlinear_channel_flow/: dummy
	+@$(def); make_lib $@

newton_solver_benchmark_set/: dummy nonlinear_channel_flow/ tosi_et_al_2015_gcubed/
	+@$(def); make_lib $@/nonlinear_channel_flow
	@$(def); run_prm $@/nonlinear_channel_flow "input_v.prm"
	@$(def); run_all_prms $@/tosi_et_al_2015
	+@$(def); make_lib $@/spiegelman_et_al_2016
	@$(def); run_prm $@/spiegelman_et_al_2016 "input.prm"


# TODO: prm doesn't run without replacing values:
onset-of-convection/: dummy
	@echo "TODO"

operator_splitting/: dummy
	+@$(def); make_lib $@/advection_reaction
	@$(def); run_all_prms $@/advection_reaction
	+@$(def); make_lib $@/exponential_decay
	@$(def); run_all_prms $@/exponential_decay

rayleigh_taylor_instability/: dummy
	+@$(def); make_lib $@
	@$(def); run_prm $@ rayleigh_taylor_instability.prm

rigid_shear/: dummy
	+@$(def); make_lib $@/plugin
	@$(def); run_all_prms $@/instantaneous
	@$(def); run_all_prms $@/time-dependent

sinking_block/: dummy
	+@$(def); make_lib $@
	@$(def); run_prm $@ sinking_block.prm

solcx/: dummy
	+@$(def); make_lib $@
	@$(def); run_all_prms $@
	+@$(def); make_lib $@/compositional_fields
	@$(def); run_all_prms $@/compositional_fields

solkz/: dummy
	+@$(def); make_lib $@
	@$(def); run_all_prms $@
	+@$(def); make_lib $@/compositional_fields
	@$(def); run_all_prms $@/compositional_fields

tangurnis/: dummy
	+@$(def); make_lib $@/code
	@$(def); run_prm $@ ba/tan.prm
	@$(def); run_prm $@ tala/tan.prm
	@$(def); run_prm $@ tala_c/tan.prm

time_dependent_annulus/: dummy
	+@$(def); make_lib $@/plugin
	@$(def); run_all_prms $@

compressibility_formulations/: dummy
	+@$(def); make_lib $@/plugins
	@$(def); run_prm $@/lateral_pipe lateral_pipe.prm
	@$(def); run_prm $@/lateral_pipe_advect lateral_pipe.prm
	@$(def); run_prm $@/lateral_pipe_increase_pressure lateral_pipe.prm
	@$(def); run_prm $@/lateral_pipe_transient lateral_pipe.prm
	@$(def); run_prm $@/vertical_pipe vertical_pipe.prm

viscoelastic_plastic_shear_bands/: dummy
	@$(def); run_all_prms $@/gerya_2019
	@$(def); run_all_prms $@/kaus_2010
