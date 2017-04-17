#!/bin/bash

# Number of processors to run on
NP=1

# averaging scheme:
for avg in "none" "arithmetic average" "geometric average" "harmonic average"
do
  echo "----Averaging scheme: $avg----"


  # global refinement:
  echo "----Global refinement----"
  for r in "4" "5" "6" "7" "8" "9"
  do
    echo "ref $r:"
    cp global.prm.base temp.prm
    echo "set Output directory = output/${avg// /_}/global/ref$r" >> temp.prm
    echo "subsection Material model" >> temp.prm
    echo "set Material averaging = $avg" >> temp.prm
    echo "end" >> temp.prm
    echo "subsection Mesh refinement" >> temp.prm
    echo "set Initial global refinement = $r" >> temp.prm
    echo "end" >> temp.prm
    mpirun -n $NP ./aspect temp.prm | grep DoFs
    rm -f temp.prm
  done


  # adaptive refinement:
  echo "----Adaptive refinement----"
  cp adaptive.prm.base temp.prm
  echo "set Output directory = output/${avg// /_}/adaptive" >> temp.prm
  echo "subsection Material model" >> temp.prm
  echo "set Material averaging = $avg" >> temp.prm
  echo "end" >> temp.prm
  mpirun -n $NP ./aspect temp.prm | grep DoFs
  rm -f temp.prm
done
