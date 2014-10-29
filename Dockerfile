FROM ubuntu:14.04

MAINTAINER Rene Gassmoeller <r.gassmoeller@mailbox.org>

ENV HOME /root

RUN apt-get update
RUN apt-get -yq install gcc \
                        build-essential \
                        wget \
                        bzip2 \
                        tar \
                        g++ \
                        gfortran \
                        libblas-dev \
                        liblapack-dev \
                        libboost-dev \
                        libopenmpi-dev \
                        openmpi-bin \
                        cmake \
                        git

#Build HDF5
RUN wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.13.tar.bz2; \
    tar xjvf hdf5-1.8.13.tar.bz2; \
    cd hdf5-1.8.13; \
    ./configure --enable-parallel \
                --enable-shared \
                --prefix=/usr/local/; \
    make -j2 && make install; \
    cd ..; \
    rm -rf /hdf5-1.8.13 /hdf5-1.8.13.tar.bz2 

#Build Trilinos
RUN wget http://trilinos.sandia.gov/download/files/trilinos-11.8.1-Source.tar.bz2; \
    tar xjvf trilinos-11.8.1-Source.tar.bz2; \
    mkdir trilinos-11.8.1-Source/build; \
    cd trilinos-11.8.1-Source/build; \
    cmake -D Trilinos_ENABLE_Sacado=ON \
          -D Trilinos_ENABLE_Stratimikos=ON \
          -D CMAKE_BUILD_TYPE=RELEASE \
          -D CMAKE_CXX_FLAGS="-O3" \
          -D CMAKE_C_FLAGS="-O3" \
          -D CMAKE_FORTRAN_FLAGS="-O5" \
          -D Trilinos_EXTRA_LINK_FLAGS="-lgfortran" \
          -D CMAKE_VERBOSE_MAKEFILE=FALSE \
          -D Trilinos_VERBOSE_CONFIGURE=FALSE \
          -D TPL_ENABLE_MPI=ON \
          -D BUILD_SHARED_LIBS=ON \
          .. ; \
    make -j4 && make install; \
    cd ../..; \
    rm -rf trilinos-11.8.1-Source trilinos-11.8.1-Source.tar.bz2 

#Build p4est
RUN wget http://p4est.github.io/release/p4est-0.3.4.2.tar.gz; \
    wget http://www.dealii.org/developer/external-libs/p4est-setup.sh; \
    chmod +x p4est-setup.sh; \
    ./p4est-setup.sh p4est-0.3.4.2.tar.gz /usr/local/p4est-0.3.4.2; \
    rm -rf p4est-build p4est-0.3.4.2 p4est-setup.sh p4est-0.3.4.2.tar.gz

#Build deal.II
RUN wget http://www.ces.clemson.edu/dealii/deal.II-8.1.0.tar.gz; \
    tar xzvf deal.II-8.1.0.tar.gz; \
    mkdir deal.II/build; \
    cd deal.II/build; \
    cmake -DDEAL_II_WITH_MPI=ON \
          -DDEAL_II_COMPONENT_COMPAT_FILES=OFF \
          -DDEAL_II_COMPONENT_EXAMPLES=OFF \
          -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DCMAKE_BUILD_TYPE=Release \
          -DP4EST_DIR=/usr/local/p4est-0.3.4.2/ \
          .. ; \
    make -j4 && make install; \
    cd ../..; \
    rm -rf deal.II deal.II-8.1.0.tar.gz

#Build aspect
RUN git clone https://github.com/geodynamics/aspect.git; \ 
    mkdir aspect/build; \
    cd aspect/build; \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DDEAL_II_DIR=/usr/local/deal.II \
          .. ; \
    make -j4; \
    mv aspect /usr/local/bin/; \
    cd /\
    rm -rf aspect

CMD ["aspect"]
