#!/bin/bash

# Number of processors to run on
NP=1
filename="output_exponential_decay"
rm $filename
echo "# Time DoFs composition_error temperature_error" >> $filename

# advection time step:
echo "----Advection time step----"
for time in "10" "5" "2.5" "1.25"
do
  echo "Advection time step $time:"
  cp exponential_decay.base.prm temp.prm
  echo "set Output directory = output/advection$time" >> temp.prm
  echo "subsection Solver parameters" >> temp.prm
  echo "subsection Operator splitting parameters" >> temp.prm
  echo "set Reaction time step = 0.0005" >> temp.prm
  echo "end" >> temp.prm
  echo "end" >> temp.prm
  echo "set Maximum time step = $time" >> temp.prm
  printf "%s ", $time | sed 's/,//g' >> $filename
  mpirun -n $NP ./aspect temp.prm | gawk '/time\=1\.000000e\+02/ {printf "%s %s %s \n", $4, $6, $10}' | sed 's/,//g' >> $filename
  rm -f temp.prm
done

# reaction time step:
echo "----Reaction time step----"
for time in "0.032" "0.016" "0.008" "0.004" "0.002" "0.001" "0.0005"
do
  echo "Reaction time step $time:"
  cp exponential_decay.base.prm temp.prm
  echo "set Output directory = output/reaction$time" >> temp.prm
  echo "subsection Solver parameters" >> temp.prm
  echo "subsection Operator splitting parameters" >> temp.prm
  echo "set Reaction time step = $time" >> temp.prm
  echo "end" >> temp.prm
  echo "end" >> temp.prm
  echo "set Maximum time step = 10" >> temp.prm
  printf "%s ", $time | sed 's/,//g' >> $filename
  mpirun -n $NP ./aspect temp.prm | gawk '/time\=1\.000000e\+02/ {printf "%s %s %s \n", $4, $6, $10}' | sed 's/,//g' >> $filename
  rm -f temp.prm
done

python ./plot_convergence.py
