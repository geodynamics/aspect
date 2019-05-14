#!/bin/bash

grid_resolution=('3' '4' '5' '6' '7' '8')
PPD=('2')

for ppd in ${PPD[@]}; do
  ppc=`expr $ppd \* $ppd`
  rm -f ppc_${ppc}.stat
  for grid_res in ${grid_resolution[@]}; do
  h=`echo "1/(2^${grid_res})" | bc -l`
  echo -n $h >> ppc_${ppc}.stat
  top_level_dir=`echo $PWD`
  current_dir=Q2_Pfalse_refinement${grid_res}
  echo $current_dir

  if [ -f ${current_dir}/log.txt ]; then
    #3 for velocity, 4 for pressure
    grep -i "Error" ${current_dir}/log.txt | cut -d':' -f2 | cut -d',' -f3 | tr -d '\n' >> ppc_${ppc}.stat
    grep -i "Error" ${current_dir}/log.txt | cut -d':' -f2 | cut -d',' -f4 >> ppc_${ppc}.stat
  fi

  done 
done
