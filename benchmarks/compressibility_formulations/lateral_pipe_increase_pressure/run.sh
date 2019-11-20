#!/bin/bash

# Run all combinations of models in this folder

repetitions=(2 4 8 16 32 64 128)
maximum_timestep=(15778800e7 7889400e7 3944700e7 1972350e7 986175e7 4930875e6 24654375e5)

for index in 0 1 2 3 4 5 6; do
  echo "subsection Geometry model" > repetitions_${repetitions[$index]}.prm
  echo "subsection Box" >> repetitions_${repetitions[$index]}.prm
  echo "set X repetitions = ${repetitions[$index]}" >> repetitions_${repetitions[$index]}.prm
  echo "end" >> repetitions_${repetitions[$index]}.prm
  echo "end" >> repetitions_${repetitions[$index]}.prm

  echo "set Maximum time step = ${maximum_timestep[$index]}" > timestep.prm
  for formulation in ala isentropic hydrostatic projected_density projected_density_full_pressure; do
    output_folder=output-lateral_pipe_increase_pressure_repetitions_${repetitions[$index]}_${formulation}_${maximum_timestep[$index]}
    echo "set Output directory = ${output_folder}" > output.prm
    cat lateral_pipe.prm ${formulation}.prm repetitions_${repetitions[$index]}.prm output.prm timestep.prm | mpirun -np 4 ../plugins/aspect --
  done
done

rm output.prm
rm repetitions_*.prm

bash make_statistics.sh
