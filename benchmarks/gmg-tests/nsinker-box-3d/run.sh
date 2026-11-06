#!/usr/bin/env bash

NP=${NP:-8}

for method in amg gmg-ls gmg-gc gmg-gc-amg; do
  for boundary_case in free-surface no-slip free-slip partial-free-slip partial-set open periodic; do
      cat "${boundary_case}.prm" > current.prm
      cat ${method}.prm.part >> current.prm

    mpirun -n "${NP}" ./aspect-release current.prm \
      | tee "log-${method}-${boundary_case}.txt"
  done
done


