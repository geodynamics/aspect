FROM ubuntu:22.04

LABEL maintainer <rene.gassmoeller@mailbox.org>

ARG VERSION=9.5.1
ARG REPO=ppa:ginggs/deal.ii-9.5.1-backports

RUN DEBIAN_FRONTEND=noninteractive apt update && apt upgrade -yq && \
  apt install -yq --install-recommends sudo software-properties-common \
  build-essential ca-certificates cmake cmake-curses-gui file gcc g++ \
  gfortran git libblas-dev liblapack-dev libopenmpi-dev \
  lsb-release ninja-build numdiff openmpi-bin \
  openmpi-common wget zlib1g-dev \
  libboost-all-dev libnetcdf-dev &&\
  add-apt-repository $REPO && apt update && \
  apt install -yq libdeal.ii-$VERSION libdeal.ii-dev && \
  apt-get clean && rm -r /var/lib/apt/lists/*

# add and enable the default user
ARG USER=dealii
ARG USER_ID=1000
RUN adduser --disabled-password --gecos '' --uid $USER_ID $USER
RUN adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#make sure everything is in place
RUN chown -R $USER:$USER /home/$USER
USER $USER
ENV HOME /home/$USER
ENV USER $USER
ENV OMPI_MCA_btl "^vader"
WORKDIR $HOME

RUN wget https://github.com/tjhei/astyle/releases/download/v2.04/astyle_2.04_linux.tar.gz && \
        tar xf astyle_2.04_linux.tar.gz && \
        cd astyle/build/gcc && make && \
        sudo make install && \
        cd && \
        rm -rf astyle* && \
        astyle --version

ENV Aspect_DIR $HOME/aspect
RUN git clone https://github.com/geodynamics/aspect.git $Aspect_DIR && \
    mkdir ${Aspect_DIR}/build && \
    cd ${Aspect_DIR}/build && \
    cmake -DCMAKE_BUILD_TYPE=DebugRelease -DCMAKE_INSTALL_PREFIX=$Aspect_DIR \
    -DASPECT_INSTALL_EXAMPLES=ON .. && \
    make -j 4 && make install && make clean && \
    ln -s $Aspect_DIR/bin/aspect $Aspect_DIR/aspect && \
    ln -s $Aspect_DIR/bin/aspect-release $Aspect_DIR/aspect-release

WORKDIR $Aspect_DIR
ENV PATH $Aspect_DIR/bin:$PATH
