#!/bin/bash
# script originated from Jonathan Perry-Houts May 2015
# Updated for installation on Stampede by D. Sarah Stamps
# at 2015 ASPECT Hackathon. Refer to stampede.pdf for
# further information.

set -x

if [ "$TRILINOS_DIR" == "" ]; then
    export TRILINOS_DIR=$HOME/packages/trilinos
    mkdir -p $TRILINOS_DIR
    echo "export TRILINOS_DIR=$TRILINOS_DIR" >> $HOME/.bashrc
fi
if [ "$P4EST_DIR" == "" ]; then
    export P4EST_DIR=$HOME/packages/p4est
    mkdir -p $P4EST_DIR
    echo "export P4EST_DIR=$P4EST_DIR" >> $HOME/.bashrc
fi
if [ "$DEAL_II_DIR" == "" ]; then
    export DEAL_II_DIR=$HOME/packages/deal.II
    mkdir -p $DEAL_II_DIR
    echo "export DEAL_II_DIR=$DEAL_II_DIR" >> $HOME/.bashrc
fi

# When installing on stampede you may need to do the following:
module load gcc/4.7.1
module load mkl/13.0.2.146
module load cmake

case $1 in

    "trilinos")
        # Build Trilinos
     if [ ! -f $HOME/Downloads/trilinos-11.12.1-Source.tar ]; then
            echo "You need trilinos-11.12.1-Source.tar in ~/Downloads"
            exit 1
        fi
        echo "Building Trilinos"
        tar -C /tmp -xf $HOME/Downloads/trilinos-11.12.1-Source.tar
        cd /tmp/trilinos-11.12.1-Source
        mkdir build ; cd build
        cmake -D Trilinos_ENABLE_Sacado=ON \
              -D Trilinos_ENABLE_Stratimikos=ON \
              -D CMAKE_BUILD_TYPE=RELEASE \
              -D CMAKE_CXX_FLAGS="-g -O3" \
              -D CMAKE_C_FLAGS="-g -O3" \
              -D CMAKE_FORTRAN_FLAGS="-g -O5" \
              -D Trilinos_EXTRA_LINK_FLAGS="-lgfortran" \
              -D CMAKE_VERBOSE_MAKEFILE=FALSE \
              -D Trilinos_VERBOSE_CONFIGURE=FALSE \
              -D TPL_ENABLE_MPI=ON \
              -D BUILD_SHARED_LIBS=ON \
              -D CMAKE_INSTALL_PREFIX=$TRILINOS_DIR \
              -D CMAKE_C_COMPILER=`which mpicc` \
              -D CMAKE_CXX_COMPILER=`which mpicxx` \
              -D CMAKE_Fortran_COMPILER=`which mpif90` \
              -D MPI_EXEC=`which mpiexec` \
              -D MPI_Fortran_COMPILER=`which mpif90` \
              -D TPL_LAPACK_LIBRARIES:STRING="/opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64" \
              -D TPL_BLAS_LIBRARIES:STRING="/opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64" \
              .. 
        make install -j7
        cd ../../ ; rm -rf trilinos-11.12.1-Source*
        ;;

    "p4est")
        # Build p4est
        echo "Building p4est"
        mkdir /tmp/p4est-build ; cd /tmp/p4est-build
        wget http://p4est.github.io/release/p4est-1.1.tar.gz
        wget http://www.dealii.org/developer/external-libs/p4est-setup.sh
        chmod +x p4est-setup.sh
        ./p4est-setup.sh p4est-1.1.tar.gz $P4EST_DIR
        cd ..; rm -rf p4est-build
        ;;

    "deal.ii")
        # Build deal.II development version
        echo "Updating deal.II"
        if [ ! -d $HOME/packages/build ]; then
           mkdir -p $HOME/packages/build
           cd $HOME/packages/build
           git clone https://github.com/dealii/dealii.git deal.II
        fi
        cd $HOME/packages/build/deal.II
        git fetch origin; git pull
        mkdir /tmp/deal.II-build ; cd /tmp/deal.II-build
        cmake -DDEAL_II_WITH_MPI=ON \
            -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
            -DTRILINOS_DIR=$TRILINOS_DIR \
            -DP4EST_DIR=$P4EST_DIR \
            -DCMAKE_CXX_FLAGS="-Wno-literal-suffix" \
            $HOME/packages/build/deal.II
        make install -j7
        cd .. ; rm -rf /tmp/deal.II-build
        ;;

    *)
        echo "Unknown option: $1"
        ;;

esac
