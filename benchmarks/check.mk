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
base_path=`pwd`; \
for file in `find . -name CMakeLists.txt`; do \
  echo "building plugin in `dirname $${file}` using ${BUILD}..."; \
  cd `dirname $${file}`; \
  rm -rf CMakeCache.txt CMakeFiles; \
  cmake -D Aspect_DIR=${BUILD} -G "Unix Makefiles" -D CMAKE_CXX_FLAGS='-Werror' . >/dev/null || { echo "cmake in `pwd` failed!"; return 1; }; \
  make >/dev/null || { echo "make in `pwd` failed!"; return 2; }; \
  echo "done building plugin in `pwd`"; \
  cd $${base_path}; \
done; \
};\
run_prm() { \
dir=$$1; \
prm=$$2; \
echo "    running $$prm in $$dir ..."; \
basepath=`pwd`; \
cd $$dir; \
[[ -f $$prm ]] || { echo "File '$$prm' missing in `pwd`."; return 7; }; \
cp $$prm $$prm.tmp || return 5; \
echo "set End time=0" >> $$prm.tmp; \
echo "set Max nonlinear iterations = 5" >> $$prm.tmp; \
${BUILD}/aspect ${CHECK} $$prm.tmp >/dev/null || { rm -f $$prm.tmp; echo "run of $$prm failed"; return 2; }; \
rm -f $$prm.tmp; \
cd $${basepath}; \
}; \
run_all_prms() { \
echo "  running all prms in $$1 ..."; \
cd $$1; \
for prm in `find . -name "*.prm" -not -path "*doc*"`; \
  do \
    run_prm `dirname $${prm}` `basename $${prm}` || return 4; \
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

# default file is too fine:
geoid-spectral-comparison/: dummy
	@$(def); run_prm $@ spectral-comparison.prm "subsection Mesh refinement\n set Initial global refinement=2\n end"

inclusion/: dummy
	+@$(def); make_lib $@
	@$(def); run_prm $@ global.prm.base
	@$(def); run_prm $@ adaptive.prm.base
	@$(def); run_all_prms $@/compositional_fields

newton_solver_benchmark_set/: dummy tosi_et_al_2015_gcubed/
	+@$(def); make_lib $@/nonlinear_channel_flow
	@$(def); run_prm $@/nonlinear_channel_flow "input_v.prm"
	@$(def); run_prm $@/nonlinear_channel_flow "input_t.prm"
	@$(def); run_all_prms $@/tosi_et_al_2015
	+@$(def); make_lib $@/spiegelman_et_al_2016
	@$(def); run_prm $@/spiegelman_et_al_2016 "input.prm"

# TODO: prm doesn't run without replacing values:
onset-of-convection/: dummy
	@echo "TODO"

tangurnis/: dummy
	+@$(def); make_lib $@/code
	@$(def); run_prm $@ ba/tan.prm
	@$(def); run_prm $@ tala/tan.prm
	@$(def); run_prm $@ tala_c/tan.prm

compressibility_formulations/: dummy
	+@$(def); make_lib $@/plugins
	@$(def); run_prm $@/lateral_pipe lateral_pipe.prm
	@$(def); run_prm $@/lateral_pipe_advect lateral_pipe.prm
	@$(def); run_prm $@/lateral_pipe_increase_pressure lateral_pipe.prm
	@$(def); run_prm $@/lateral_pipe_transient lateral_pipe.prm
	@$(def); run_prm $@/vertical_pipe vertical_pipe.prm

# amg.prm does not converge with a coarse mesh
nsinker_spherical_shell/: dummy
	+@$(def); make_lib $@
	@$(def); run_prm $@ gmg.prm
