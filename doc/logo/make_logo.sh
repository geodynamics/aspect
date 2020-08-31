#!/bin/bash

mkdir -p logo
cd logo && mpirun -np 6 ../logo-build/aspect logo.prm
