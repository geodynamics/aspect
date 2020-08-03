#!/bin/perl
#
# This script removes references to the isothermal compression formulation
# and replaces them with the isentropic compression formulation.

while (<>)
{
  s/isothermal_compression/isentropic_compression/g;
  s/isothermal compression/isentropic compression/g;
  s/IsothermalCompression/IsentropicCompression/g;
  print;
}
