#!/bin/bash

for directory in `ls . | grep pipe`; do
  cd $directory
  bash run.sh
  cd ..
done
