#!/bin/bash

refinements=('2' '3' '4')
thermal_conductivities=('1e-2' '1e-3')
temperature_degrees=('1' '2' '3')

for temperature_degree in ${temperature_degrees[@]}; do
  for thermal_conductivity in ${thermal_conductivities[@]}; do
    outputfile=results_Q${temperature_degree}_${thermal_conductivity}.dat
    rm -f $outputfile

      for refinement in ${refinements[@]}; do
        h=`echo "1/(2^${refinement})" | bc -l`

        current_dir=Q${temperature_degree}_refinement${refinement}_${thermal_conductivity}
        echo $current_dir
        echo -n $h >> $outputfile
        echo -n ' ' >> $outputfile

        if [ -f ${current_dir}/log.txt ]; then
          # 1 average temperature, 2 maximum temperature
          grep -i "Temperature min/avg/max" ${current_dir}/log.txt | tail -n 1 | cut -d':' -f2 | cut -d',' -f2 | cut -d' ' -f2 | tr -d '\n' >> $outputfile
          echo -n  ' ' >> $outputfile
          grep -i "Temperature min/avg/max" ${current_dir}/log.txt | tail -n 1 | cut -d':' -f2 | cut -d',' -f3 | cut -d' ' -f2 | tr -d '\n' >> $outputfile
        fi
	echo '' >> $outputfile

      done
  done 
done
