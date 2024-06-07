#!/bin/bash

output_file="timings.csv"
echo "refinement, model_time_euler, particle_advect_euler, advect_percent_euler," \
    "model_time_rk2, particle_advect_rk2, advect_percent_rk2," \
    "model_time_rk4, particle_advect_rk4, advect_percent_rk4," \
    "model_time_rk4_38, particle_advect_rk4_38, advect_percent_rk4_38," \
    "rk2_overhead" > $output_file

for refinement in 2 3 4 5; do # 
  # velocity_accuracy_colum=`grep "u_L2" refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/statistics | cut -d ' ' -f 2 | head -c -2`
  # velocity_error_euler=`tail -n 1 refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/statistics | cut -d ' ' -f $velocity_accuracy_colum`
  # velocity_error_rk2=`tail -n 1 refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/statistics | cut -d ' ' -f $velocity_accuracy_colum`
  # velocity_error_rk4=`tail -n 1 refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/statistics | cut -d ' ' -f $velocity_accuracy_colum`
  # velocity_error_rk4_38=`tail -n 1 refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/statistics | cut -d ' ' -f $velocity_accuracy_colum`



  # replace exponential notation with 10^
  # velocity_error_euler=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$velocity_error_rk2"`
  # velocity_error_rk2=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$velocity_error_rk2"`
  # velocity_error_rk4=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$velocity_error_rk2"`
  # velocity_error_rk4_38=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$velocity_error_rk2"`


  model_time_euler=`grep "Total wallclock time elapsed since start" refinement_${refinement}_higher_order_false_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 3 | head -c -3`
  model_time_rk2=`grep "Total wallclock time elapsed since start" refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 3 | head -c -3`
  model_time_rk4=`grep "Total wallclock time elapsed since start" refinement_${refinement}_higher_order_false_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 3 | head -c -3`
  model_time_rk4_38=`grep "Total wallclock time elapsed since start" refinement_${refinement}_higher_order_false_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 3 | head -c -3`


  # replace exponential notation with 10^
  #model_time_euler=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$model_time_euler"`  
  #model_time_rk2=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$model_time_rk2"`
  #model_time_rk4=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$model_time_rk4"`
  #model_time_rk4_38=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$model_time_rk4_38"`


  particle_time_euler=`grep "Particles: Advect" refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 4 | head -c -3`
  particle_time_rk2=`grep "Particles: Advect" refinement_${refinement}_higher_order_true_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 4 | head -c -3`
  particle_time_rk4=`grep "Particles: Advect" refinement_${refinement}_higher_order_false_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 4 | head -c -3`
  particle_time_rk4_38=`grep "Particles: Advect" refinement_${refinement}_higher_order_false_interpolation_bilinear_least_squares/log.txt | cut -d \| -f 4 | head -c -3`


  particle_percentage_euler=`echo "scale=3; $particle_time_euler / $model_time_rk2 * 100" | bc`
  particle_percentage_rk4=`echo "scale=3; $particle_time_rk4 / $model_time_rk2 * 100" | bc`
  particle_percentage_rk4_38=`echo "scale=3; $particle_time_rk4_38 / $model_time_rk2 * 100" | bc`


  euler_overhead=`echo "scale=3; ($particle_time_euler / $particle_time_rk2 - 1.0) * 100" | bc`
  rk4_overhead=`echo "scale=3; ($particle_time_rk4 / $particle_time_rk2 - 1.0) * 100" | bc`
  rk4_38_overhead=`echo "scale=3; ($particle_time_rk4_38 / $particle_time_rk2 - 1.0) * 100" | bc`


  # euler_accuracy_per_time=`echo "scale=10; ((1.0-($velocity_error_euler)) / $model_time_rk2)" | bc`
  # rk2_accuracy_per_time=`echo "scale=10; (1.0-($velocity_error_rk2)) / $model_time_rk2" | bc`
  # rk4_accuracy_per_time=`echo "scale=10; ((1.0-($velocity_error_rk4)) / $model_time_rk2)" | bc`
  # rk4_38_accuracy_per_time=`echo "scale=10; ((1.0-($velocity_error_rk4_38)) / $model_time_rk2" | bc`


  echo "$refinement, $model_time_rk2, $particle_time_rk2, $particle_percentage_rk2," \
      "$model_time_rk4, $particle_time_rk4, $particle_percentage_rk4, $rk2_overhead" >> $output_file
done


