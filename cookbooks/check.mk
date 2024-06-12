#
# This Makefile compiles all plugins and runs all prm files in the
# subdirectories of the cookbooks folder
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
for file in *.prm; \
  do \
    if [[ "$${file}" != *.part.prm ]]; then \
      run_prm . $${file} || return 4; \
    fi \
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

mainprms:= $(wildcard *.prm)

%.prm: dummy
	@$(def); run_prm . $@

main: dummy $(mainprms)


#
#
# Custom rules. Make them dependent on dummy like this:
#
# example/: dummy
#	@$(def); run_prm $@ test.prm

# This does not run without generating an input file using a python script
prescribed_velocity_ascii_data/: dummy
	+@$(def); make_lib $@

# Fastscape is not compiled by default, so we cannot test this
fastscape_eroding_box/: dummy
	

# Manually list files (and ignore mantle_setup_restart.prm as _start does not generate
# a checkpoint file after 1 step
future/: dummy
	@$(def); run_prm $@ mantle_setup_start.prm
	@$(def); run_prm $@ net_rotation.prm
	@$(def); run_prm $@ periodic_box.prm
	@$(def); run_prm $@ radiogenic_heating.prm
	@$(def); run_prm $@ radiogenic_heating_function.prm
	@$(def); run_prm $@ sphere.prm
