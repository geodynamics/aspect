#!/bin/bash
#
# This script compiles all plugins and run all prm files in the subdirectories
# of the benchmark folder

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


# run aspect on all .prm files in the current folder or any subdirectory
run_prm ()
{
    for prm in `find . -name "*prm"`;
    do
    if [ "`basename $prm`" = "parameters.prm" ];
    then
	continue;
    fi
    echo "Running '$prm' at `pwd` with '$BUILD' ..."
    cp $prm tmp.prm~
    echo "set End time=0" >> tmp.prm~
    $BUILD/aspect tmp.prm~ >/dev/null || { rm -f tmp.prm~; return 2; }
    rm -f tmp.prm~
    done
    return 0;
}

# configure and compile the plugin in the current directory
make_lib ()
{
    echo "configuring in `pwd` ..."
    rm -rf CMakeCache.txt
    cmake -D ASPECT_DIR=$BUILD . >/dev/null || { echo "cmake failed!"; return 1; }
    make >/dev/null || { echo "make failed!"; return 2; }
    return 0;
}


echo "Checking benchmarks using $BUILD/aspect ..."

(cd burstedde; make_lib && run_prm ) || { echo "FAILED"; exit 1; } 

(cd davies_et_al; cd case-2.3-plugin; make_lib && cd .. && run_prm ) || { echo "FAILED"; exit 1; } 

(cd inclusion; make_lib && run_prm ) || { echo "FAILED"; exit 1; } 

(cd solcx; make_lib && run_prm ) || { echo "FAILED"; exit 1; } 

(cd solkz; make_lib && run_prm ) || { echo "FAILED"; exit 1; } 

(cd tangurnis; cd code; make_lib && cd .. && run_prm ) || { echo "FAILED"; exit 1; } 

echo "all good! :-)"
exit 0
