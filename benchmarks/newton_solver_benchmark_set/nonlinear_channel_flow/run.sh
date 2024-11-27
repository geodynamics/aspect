#!/bin/bash

processes=4
outfile="output.log"
errorfile="error.log"
plotfile="plot.dat"

for boundary_type in "v" "t"; do
  infile="input_${boundary_type}.prm"

  # First run one model with iterated Advection and Stokes
  dirname_clean="BT_${boundary_type}_singleAdvectioniteratedStokes"
  dirname="results/${dirname_clean}"
  mkdir -p $dirname
  echo "Starting $dirname"

  sed  \
    -e "s/set Output directory .*/set Output directory = results\/$dirname_clean/g" \
    -e "s/set Nonlinear solver scheme.*/set Nonlinear solver scheme = single Advection, iterated Stokes/g" \
    -e "s/set List of output variables .*/set List of output variables = material properties,strain rate/g" \
    -e "s/set Linear solver tolerance .*/set Linear solver tolerance = 1e-8/g" \
    $infile > $dirname/$infile

  nohup mpirun -np $processes ./aspect $dirname/$infile > $dirname/$outfile 2>$dirname/$errorfile
  grep "Relative nonlinear residual " $dirname/$outfile > $dirname/$plotfile

  # Then run all the Newton solver models
  for stabilization in "SPD" "none"; do
    for picard_iterations in 5 10 15 150; do
      for line_search_iterations in 0 5 100; do
        for max_linear_tolerance in "1e-8" "1e-2"; do
          for residual_scaling_method in "false" "true"; do
            dirname_clean="BT_${boundary_type}_singleAdvectioniteratedNewtonStokes_stabilization_${stabilization}_mLT_${max_linear_tolerance}_P_${picard_iterations}_LS_${line_search_iterations}_RSM_${residual_scaling_method}"
            dirname="results/${dirname_clean}"
            mkdir -p $dirname
            echo "Starting $dirname"

            sed  \
                -e "s/set Output directory .*/set Output directory = results\/$dirname_clean/g" \
                -e "s/set Max pre-Newton nonlinear iterations = .*/set Max pre-Newton nonlinear iterations = $picard_iterations/g" \
                -e "s/set Stabilization preconditioner = .*/set Stabilization preconditioner = $stabilization/g" \
                -e "s/set Stabilization velocity block = .*/set Stabilization velocity block = $stabilization/g" \
                -e "s/set Max Newton line search iterations = .*/set Max Newton line search iterations = $line_search_iterations/g" \
                -e "s/set Maximum linear Stokes solver tolerance =.*/set Maximum linear Stokes solver tolerance = $max_linear_tolerance/g" \
                -e "s/set Use Newton residual scaling method .*/set Use Newton residual scaling method = $residual_scaling_method/g" \
                $infile > $dirname/$infile

            nohup mpirun -np $processes ./aspect $dirname/$infile > $dirname/$outfile 2>$dirname/$errorfile
            grep "Relative nonlinear residual " $dirname/$outfile > $dirname/$plotfile
          done
        done
      done
    done
  done
done
