#!/bin/bash

# install directory (if relative, relative to current working directory)
install_dir="./"

install_dir=`python -c "import os.path; print os.path.abspath('${install_dir}')"`

# download perplex sources from the internet if necessary
if [ ! -f ./perplex-sources-stable.zip ]; then
    wget https://petrol.natur.cuni.cz/~ondro/perplex-sources-stable.zip
fi

# unpack zip file (if necessary)
if [ ! -d ./perplex-sources ]; then
    unzip perplex-sources-stable.zip
fi

# Remove zip file now that everything is unpacked
rm perplex-sources-stable.zip

# Change directory name
mv perplex-stable source

# unpack cmake build
cd source
cp ../cmake_perplex.tar .
tar xvf cmake_perplex.tar

# make c header
cd include
python ./create_c_header.py
cd ..

# set up cmake, install to chosen install directory
cmake -D CMAKE_INSTALL_PREFIX=${install_dir}/install .

# install perplex to chosen directory
make -fMakefile install

cd ..


