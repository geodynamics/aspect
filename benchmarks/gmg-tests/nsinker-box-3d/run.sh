#!/usr/bin/env bash

NP=${NP:-8}

for method in GMG-GC GMG-LS AMG; do
  for boundary_case in free-surface no-slip free-slip partial-free-slip partial-set open periodic; do
    cat "${boundary_case}.prm" > current.prm

    if [[ ${method} == AMG ]]; then
      cat amg.prm.part >> current.prm
    elif [[ ${method} == GMG-LS ]]; then
      cat gmg-local-smoothing.prm.part >> current.prm
    fi

    mpirun -n "${NP}" ./aspect-release current.prm \
      | tee "log-${method}-${boundary_case}.txt"
  done
done


