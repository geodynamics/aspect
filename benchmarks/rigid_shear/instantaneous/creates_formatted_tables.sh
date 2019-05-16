#!/bin/bash

grid_resolution=('2' '3' '4' '5' '6' '7' '8')
stokes_degrees=('2' '3')

for stokes_degree in ${stokes_degrees[@]}; do
  rm -f results.txt
  for grid_res in ${grid_resolution[@]}; do
    h=`echo "1/(2^${grid_res})" | bc -l`
    echo -n $h >> results.txt

    for stokes_degree in ${stokes_degrees[@]}; do
      current_dir=Q${stokes_degree}_Pfalse_refinement${grid_res}
      echo $current_dir

      if [ -f ${current_dir}/log.txt ]; then
        #3 for velocity, 4 for pressure, 5 for density
        grep -i "Error" ${current_dir}/log.txt | cut -d':' -f2 | cut -d',' -f3 | tr -d '\n' >> results.txt
        grep -i "Error" ${current_dir}/log.txt | cut -d':' -f2 | cut -d',' -f4 | tr -d '\n' >> results.txt
        grep -i "Error" ${current_dir}/log.txt | cut -d':' -f2 | cut -d',' -f5 | tr -d '\n' >> results.txt
      fi
    done
    echo "" >> results.txt
  done 
done
