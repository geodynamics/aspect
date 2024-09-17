#!/bin/bash

for int in euler rk2 rk4; do #euler rk2 rk4
  rm max_deviation_$int

  for cfl in 1.0 0.5 0.25 0.125 0.0625 0.03125; do
  # Prepare and run model
  OUTPUT_DIR=circle_${int}_${cfl}
  echo set Output directory = $OUTPUT_DIR > output_dir.prm
  echo set CFL number = $cfl > cfl.prm
  
  cat circle.prm output_dir.prm cfl.prm ${int}.part.prm | aspect --

  # format and print results
    tail -n +8 $OUTPUT_DIR/particles/particles-00000.0000.gnuplot | sort -g -k 3 > sort_temp_0
    tail -n +8 $OUTPUT_DIR/particles/particles-00001.0000.gnuplot | sort -g -k 3 > sort_temp_1
    paste sort_temp_0 sort_temp_1 > sort_temp
    gawk -v cfl=$cfl 'BEGIN{max_deviation=0.0} {deviation=sqrt(($1-$4)*($1-$4)+($2-$5)*($2-$5));if(deviation>max_deviation) {max_deviation=deviation}} END{print cfl, max_deviation}' sort_temp >> max_deviation_$int

  
  done
done

# make visualization
echo set Output directory = circle_vis > output_dir.prm
cat circle.prm vis.prm output_dir.prm | aspect --

# make plot
#gnuplot plot_convergence.plt

rm output_dir.prm cfl.prm sort_temp_0 sort_temp_1 sort_temp
