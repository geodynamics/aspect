#!/bin/bash

# Number of processors to run on
NP=1
filename="output_reaction_advection_high_res"
rm $filename
echo "# Time DoFs composition_error temperature_error" >> $filename

# advection time step:
echo "----Advection time step----"
for time in "1" "0.5" "0.25" "0.125" "0.0625" "0.03125" "0.015625"
do
  echo "Advection time step $time:"
  cp advection_reaction.base.prm temp.prm
  echo "set Output directory = output/advection$time" >> temp.prm
  echo "subsection Solver parameters" >> temp.prm
  echo "subsection Operator splitting parameters" >> temp.prm
  echo "set Reaction time step = 0.001953125" >> temp.prm
  echo "end" >> temp.prm
  echo "end" >> temp.prm
  echo "set Maximum time step = $time" >> temp.prm
  printf "%s ", $time | sed 's/,//g' >> $filename
  mpirun -n $NP ./aspect temp.prm | gawk '/time\=1\.000000e\+00/ {printf "%s %s %s \n", $4, $6, $10}' | sed 's/,//g' >> $filename
  rm -f temp.prm
done

# reaction time step:
echo "----Reaction time step----"
for time in "0.0125" "0.0064" "0.0032" "0.0016" "0.0008"
do
  factor="10"
  echo "Reaction time step $time:"
  cp advection_reaction.base.prm temp.prm
  echo "set Output directory = output/reaction$time" >> temp.prm
  echo "subsection Solver parameters" >> temp.prm
  echo "subsection Operator splitting parameters" >> temp.prm
  echo "set Reaction time step = $time" >> temp.prm
  echo "end" >> temp.prm
  echo "end" >> temp.prm
  advection_time="$(echo $time*$factor | bc)"
  echo "set Maximum time step = $advection_time" >> temp.prm
  printf "%s ", $time | sed 's/,//g' >> $filename
  mpirun -n $NP ./aspect temp.prm | gawk '/time\=1\.000000e\+00/ {printf "%s %s %s \n", $4, $6, $10}' | sed 's/,//g' >> $filename
  rm -f temp.prm
done

python ./plot_convergence_high_res.py
