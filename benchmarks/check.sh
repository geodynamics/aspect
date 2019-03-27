#!/bin/bash
#
# This script compiles all plugins and runs all prm files in the subdirectories
# of the benchmarks folder

if [ "$#" -ne 1 ]; then
    echo "usage: $0 aspect-build-directory"
    exit 1
fi

BUILD=`cd $1;pwd`

if [[ ! -e $BUILD/AspectConfig.cmake || ! -e $BUILD/aspect  ]];
then
    echo "'$BUILD' doesn't look like a build directory"
    exit 1
fi

#
run_prm ()
{
    for prm in "$@"
    do
        echo "Running '$prm' at `pwd` with '$BUILD' ..."
        cp $prm $prm.tmp
        echo "set End time=0" >> $prm.tmp
        $BUILD/aspect $prm.tmp >/dev/null || { rm -f $prm.tmp; return 2; }
        rm -f $prm.tmp
    done
}

# run aspect on all .prm files in the current folder or any subdirectory
run_all_prms ()
{
    for prm in `find . -name "*prm"`;
    do
        if [ "`basename $prm`" = "parameters.prm" ];
        then
	        continue;
        fi
        echo "Running '$prm' at `pwd` with '$BUILD' ..."
        cp $prm $prm.tmp
        echo "set End time=0" >> $prm.tmp
        $BUILD/aspect $prm.tmp >/dev/null || { rm -f $prm.tmp; return 2; }
        rm -f $prm.tmp
    done
    echo "... completed `pwd`"
    return 0;
}

# configure and compile the plugin in the current directory
make_lib ()
{
    echo "configuring in `pwd` ..."
    rm -rf CMakeCache.txt
    cmake -D Aspect_DIR=$BUILD . >/dev/null || { echo "cmake failed!"; return 1; }
    make >/dev/null || { echo "make failed!"; return 2; }
    return 0;
}


echo "Checking benchmarks using $BUILD/aspect"
echo "Please be patient..."

if false; then

( (cd annulus; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd blankenbach/plugin; make_lib && cd .. && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

#TODO: broken
#( (cd buiter_et_al_2008_jgr; run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd burstedde; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd crameri_et_al/case_1 && make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd crameri_et_al/case_2 && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd davies_et_al; cd case-2.3-plugin; make_lib && cd .. && run_prm "case-2.1.prm" && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd doneahuerta/; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd finite_strain && make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd geoid-spectral-comparison; run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd hollow_sphere; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd inclusion; make_lib && run_prm "global.prm.base" "adaptive.prm.base") || { echo "FAILED"; exit 1; } ) &

wait

( (cd inclusion/compositional_fields; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd king2dcompressible && make_lib && run_prm "ala.prm" ) || { echo "FAILED"; exit 1; } ) &

( (cd layeredflow && make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

# TODO: no .prm in folder:
( (cd nonlinear_channel_flow && make_lib ) || { echo "FAILED"; exit 1; } ) &

wait # newton_solver_benchmark_set/nonlinear_channel_flow depends on nonlinear_channel_flow/

( (cd newton_solver_benchmark_set/nonlinear_channel_flow/ && run_prm "input_v.prm" ) || { echo "FAILED"; exit 1; } ) &

# TODO: broken
#( (cd newton_solver_benchmark_set/tosi_et_al_2015/ && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd nsinker && make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

# TODO: prm doesn't run without replacing values:
#( (cd onset-of-convection && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd operator_splitting/advection_reaction && make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd operator_splitting/exponential_decay/ && make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd shear_bands; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd solcx; make_lib && run_prm "solcx.prm" ) || { echo "FAILED"; exit 1; } ) &

( (cd solcx/compositional_fields; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

fi

( (cd solitary_wave; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd solkz; make_lib && run_prm "solkz.prm" ) || { echo "FAILED"; exit 1; } ) &

( (cd solkz/compositional_fields; make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd tangurnis; cd code; make_lib && cd .. && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd time_dependent_annulus/plugin && make_lib && cd .. && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd tosi_et_al_2015_gcubed/ && make_lib && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd viscoelastic_bending_beam && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd viscoelastic_stress_build-up && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

( (cd zhong_et_al_93 && run_all_prms ) || { echo "FAILED"; exit 1; } ) &

wait

exit 0
