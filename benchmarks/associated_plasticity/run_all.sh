#!/bin/bash

mpirun -n 4 ./aspect implicit.prm
mpirun -n 4 ./aspect explicit.prm
mpirun -n 4 ./aspect disabled.prm
