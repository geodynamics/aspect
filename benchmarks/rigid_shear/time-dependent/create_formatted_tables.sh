#!/bin/bash

grid_resolution=('3' '4' '5' '6')
interpolators=('cell_average' 'bilinear_least_squares')
integrators=('rk2' 'rk4')
stokes_degrees=('2' '3')

for stokes_degree in ${stokes_degrees[@]}; do
  for interpolator in ${interpolators[@]}; do
    for integrator in ${integrators[@]}; do
      outputfile=results_Q${stokes_degree}_${interpolator}_${integrator}.dat
      rm -f $outputfile
      echo "h                     ppc e_u_L1      e_p_L1       e_u_L2       e_p_L2       e_rho_L2" > $outputfile

      if [ $stokes_degree == 2 ]; then
        PPD=('4' '5' '6' '7' '10' '15' '20')
      else
        PPD=('4' '5' '6' '7' '10' '15' '20' '32' '45' '64' '80')
      fi

      for ppd in ${PPD[@]}; do
        for grid_res in ${grid_resolution[@]}; do
          h=`echo "1/(2^${grid_res})" | bc -l`
          ppc=`expr $ppd \* $ppd`
          echo -n $h $ppc >> $outputfile

          current_dir=Q${stokes_degree}_Pfalse_refinement${grid_res}_${ppd}_reference_cell_${interpolator}_${integrator}
          echo $current_dir

          if [ -f ${current_dir}/log.txt ]; then
            # 1 L1_velocity, 2 L1_pressure, 3 L2_velocity, 4 L2_pressure, 5 L2_density
            grep -i "Error" ${current_dir}/log.txt | tail -n 1 | cut -d':' -f2 | cut -d',' -f1 | tr -d '\n' >> $outputfile
            grep -i "Error" ${current_dir}/log.txt | tail -n 1 | cut -d':' -f2 | cut -d',' -f2 | tr -d '\n' >> $outputfile
            grep -i "Error" ${current_dir}/log.txt | tail -n 1 | cut -d':' -f2 | cut -d',' -f3 | tr -d '\n' >> $outputfile
            grep -i "Error" ${current_dir}/log.txt | tail -n 1 | cut -d':' -f2 | cut -d',' -f4 | tr -d '\n' >> $outputfile
            grep -i "Error" ${current_dir}/log.txt | tail -n 1 | cut -d':' -f2 | cut -d',' -f5 | tr -d '\n' >> $outputfile
          fi
	  echo '' >> $outputfile

        done
      done
    done
  done 
done
