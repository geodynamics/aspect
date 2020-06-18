#!/bin/bash

git clone https://github.com/geodynamics/aspect.git

cd aspect
git checkout 148bb56
cd ..

mkdir logo-build
cd logo-build
cmake ../aspect
make release
make -j6
cd ..
