#!/bin/perl
#
# This script uses the new n_evaluation_points functions
# of the material model inputs and outputs objects

while (<>)
{
  s/in.position.size/in.n_evaluation_points/g;
  s/in.temperature.size/in.n_evaluation_points/g;
  s/out.viscosities.size/out.n_evaluation_points/g;
  print;
}
