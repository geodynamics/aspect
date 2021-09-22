#!/bin/perl
#
# This script removes remnants of the parameters.prm
# output, in particular the default value appended to lines.

while (<>)
{
  s/\s*# default:.*$//g;
  print;
}
