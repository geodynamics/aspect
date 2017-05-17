#!/bin/bash

# Number of processors to run on
NP=1


# advection time step:
echo "----Advection time step----"
for time in "10" "5" "2.5" "1.25"
do
  echo "Advection time step $time:"
  cp advection_reaction.base.prm temp.prm
  echo "set Output directory = output/advection$time" >> temp.prm
  echo "set Reaction time step = 0.0005" >> temp.prm
  echo "set Maximum time step = $time" >> temp.prm
  mpirun -n $NP ./aspect temp.prm | grep "time=1.000000e+02"
  rm -f temp.prm
done

# reaction time step:
echo "----Reaction time step----"
for time in "0.032" "0.016" "0.008" "0.004" "0.002" "0.001" "0.0005"
do
  echo "Reaction time step $time:"
  cp advection_reaction.base.prm temp.prm
  echo "set Output directory = output/reaction$time" >> temp.prm
  echo "set Reaction time step = $time" >> temp.prm
  echo "set Maximum time step = 10" >> temp.prm
  mpirun -n $NP ./aspect temp.prm | grep "time=1.000000e+02"
  rm -f temp.prm
done
