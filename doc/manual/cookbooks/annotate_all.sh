#!/bin/bash

# this script applies annotate.pl to all *prm in the current and all
# subdirectories and creates *prm.out files that are then used in the manual.

find . -name \*.prm|while read fname; do
  echo $fname
  perl annotate.pl "$fname" >"$fname.out"
done
