FROM ubuntu:14.04

LABEL maintainer <rene.gassmoeller@mailbox.org>

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

RUN adduser aspect_user
USER aspect_user
ENV HOME /home/aspect_user
WORKDIR /home/aspect_user

#Build HDF5
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.0-patch1/src/hdf5-1.10.0-patch1.tar.bz2; \
    tar xjvf hdf5-1.10.0-patch1.tar.bz2; \
    cd hdf5-1.10.0-patch1; \
    ./configure --enable-parallel \
                --enable-shared \
                --prefix=$HOME/hdf5-1.10.0-install; \
    make -j4 && make install; \
    cd ..; \
    rm -rf hdf5-1.10.0-patch1 hdf5-1.10.0-patch1.tar.bz2
ENV HDF5_DIR /home/aspect_user/hdf5-1.10.0-install

#Build Trilinos
RUN wget http://trilinos.sandia.gov/download/files/trilinos-11.8.1-Source.tar.bz2; \
    tar xjvf trilinos-11.8.1-Source.tar.bz2; \
    mkdir trilinos-11.8.1-Source/build; \
    cd trilinos-11.8.1-Source/build; \
    cmake -D Trilinos_ENABLE_Sacado=ON \
          -D Trilinos_ENABLE_Stratimikos=ON \
          -D CMAKE_BUILD_TYPE=RELEASE \
          -D CMAKE_INSTALL_PREFIX=$HOME/trilinos-11.8.1-install \
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
ENV TRILINOS_DIR /home/aspect_user/trilinos-11.8.1-install

#Build p4est
RUN wget http://p4est.github.io/release/p4est-0.3.4.2.tar.gz; \
    wget http://www.dealii.org/developer/external-libs/p4est-setup.sh; \
    chmod +x p4est-setup.sh; \
    ./p4est-setup.sh p4est-0.3.4.2.tar.gz $HOME/p4est-0.3.4.2-install; \
    rm -rf p4est-build p4est-0.3.4.2 p4est-setup.sh p4est-0.3.4.2.tar.gz
ENV P4EST_DIR /home/aspect_user/p4est-0.3.4.2-install

#Build deal.II
RUN wget http://github.com/dealii/dealii/releases/download/v8.4.2/dealii-8.4.2.tar.gz && \
    tar xzvf dealii-8.4.2.tar.gz && \
    mkdir dealii-8.4.2/build && \
    cd dealii-8.4.2/build && \
    cmake -DDEAL_II_WITH_MPI=ON \
          -DDEAL_II_COMPONENT_COMPAT_FILES=OFF \
          -DDEAL_II_COMPONENT_EXAMPLES=OFF \
          -DCMAKE_INSTALL_PREFIX=$HOME/deal.II-8.4.2-install \
          -DCMAKE_BUILD_TYPE=Release \
          .. && \
    make -j4 && make install && \
    cd ../.. && \
    rm -rf dealii-8.4.2 dealii-8.4.2.tar.gz
ENV DEAL_II_DIR /home/aspect_user/deal.II-8.4.2-install

#Build aspect
RUN git clone https://github.com/geodynamics/aspect.git ./aspect && \ 
    mkdir aspect/build && \
    cd aspect/build && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DDEAL_II_DIR=$HOME/deal.II-install \
          .. && \
    make -j4 && \
    mv aspect $HOME/aspect/ && \
    cd .. && \
    rm -rf build

WORKDIR /home/aspect_user/aspect
