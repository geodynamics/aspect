#!/bin/bash

find . -name \*.prm|while read fname; do
  echo $fname
  perl annotate.pl "$fname" >"$fname.out"
done
