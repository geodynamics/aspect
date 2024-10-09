#!/bin/perl
#
# This script removes references to the particle world class and replaces
# them with the new particle manager class.

while (<>)
{
  s/get_particle_world\(\)/get_particle_manager\(0\)/g;

  print;
}
