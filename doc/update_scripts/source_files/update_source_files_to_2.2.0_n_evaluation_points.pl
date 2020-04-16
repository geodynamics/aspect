#!/bin/perl
#
# This script removes uses the new n_evaluation_points functions
# of the material model inputs object

while (<>)
{
  s/in.position.size/in.n_evaluation_points/g;
  s/in.temperature.size/in.n_evaluation_points/g;
  print;
}
