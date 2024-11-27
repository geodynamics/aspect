#!/bin/bash
processes=4
infile="input.prm"
outfile="output.log"
errorfile="error.log"
plotfile="plot.dat"

for case in 1 2 3; do
  # velocity in m/s
  if [ $case == 1 ]; then
    velocity=25
    U="7.92219116e-11"
    background_viscosity="1e23"
  elif [ $case == 2 ]; then
    velocity=50
    U="1.58443823e-10"
    background_viscosity="1e24"
  elif [ $case == 3 ]; then
    velocity=125
    U="3.96109558e-10"
    background_viscosity="5e24"
  fi

  # First run one model with only Picard iterations
  dirname_base="singleAdvectionanditeratedStokes_vel_${velocity}_BV_${background_viscosity}"
  dirname="results/${dirname_base}"
  echo "Starting $dirname"
  mkdir -p $dirname

  sed  \
    -e "s/set Function expression = if(x<60e3,.*/set Function expression = if(x<60e3,$U,-$U);0/g" \
    -e "s/set Reference viscosity .*/set Reference viscosity = $background_viscosity/g" \
    -e "s/set Output directory .*/set Output directory = results\/$dirname_base/g" \
    -e "s/set Nonlinear solver scheme .*/set Nonlinear solver scheme = single Advection, iterated Stokes/g" \
    -e "s/set List of output variables .*/set List of output variables = material properties,strain rate/g" \
    -e "s/set Linear solver tolerance .*/set Linear solver tolerance = 1e-14/g" \
    $infile > "$dirname/$infile"

  nohup mpirun -np $processes ./aspect "${dirname}/${infile}" > ${dirname}/${outfile} 2>${dirname}/${errorfile}
  grep "Relative nonlinear residual " ${dirname}/${outfile} > ${dirname}/${plotfile}

  # Now run Newton solver models with varying solver parameters
  for stabilization in "none" "SPD"; do
    for picard_iterations in 0 5 15 25 150; do # note 150 represents pure defect-correction Picard
      for maximum_linear_tolerance in "9e-1" "1e-8"; do
        dirname_base="singleAdvectionandNewtonStokes""_NS_${stabilization}_mLT_${maximum_linear_tolerance}_P${picard_iterations}_vel_${velocity}_BV_${background_viscosity}"
        dirname="results/$dirname_base"
        echo "Starting $dirname"
        mkdir -p $dirname

        sed  \
          -e "s/set Function expression = if(x<60e3,.*/set Function expression = if(x<60e3,$U,-$U);0/g" \
          -e "s/set Reference viscosity .*/set Reference viscosity = $background_viscosity/g" \
          -e "s/set Output directory .*/set Output directory = results\/$dirname_base/g" \
          -e "s/set Max pre-Newton nonlinear iterations .*/set Max pre-Newton nonlinear iterations = $picard_iterations/g" \
          -e "s/set Stabilization preconditioner .*/set Stabilization preconditioner = $stabilization/g" \
          -e "s/set Stabilization velocity block .*/set Stabilization velocity block = $stabilization/g" \
          -e "s/set Maximum linear Stokes solver tolerance .*/set Maximum linear Stokes solver tolerance = $maximum_linear_tolerance/g" \
          $infile > "$dirname/$infile"

        nohup mpirun -np $processes ./aspect $dirname/$infile > $dirname/$outfile 2>$dirname/$errorfile
        grep "Relative nonlinear residual " $dirname/$outfile > $dirname/$plotfile
      done
    done
  done
done
