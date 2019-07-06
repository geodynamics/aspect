#!/bin/bash

# Run all combinations of models in this folder

for repetitions in 2 4 8 16 32 64 128; do #1 2 4 8 16 32 64 128
  echo "subsection Geometry model" > repetitions.prm
  echo "subsection Box" >> repetitions.prm
  echo "set X repetitions = $repetitions" >> repetitions.prm
  echo "end" >> repetitions.prm
  echo "end" >> repetitions.prm

  for temperature in sub_adiabatic adiabatic; do
    for formulation in ala isentropic hydrostatic projected_density; do
      output_folder=output-lateral_pipe_repetitions_${repetitions}_${formulation}_${temperature}
      echo "set Output directory = ${output_folder}" > output.prm
      cat lateral_pipe.prm ${temperature}.prm ${formulation}.prm repetitions.prm output.prm | mpirun -np 4 ../plugins/aspect --
    done
  done
done

rm output.prm
rm repetitions.prm

bash make_statistics.sh
