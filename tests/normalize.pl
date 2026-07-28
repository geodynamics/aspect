#!/usr/bin/perl
######################################################################
#
# Copyright (C) 2017 - 2026 by the authors of the ASPECT code.
#
#  This file is part of ASPECT.
#
#  ASPECT is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2, or (at your option)
#  any later version.
#
#  ASPECT is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ASPECT; see the file LICENSE.  If not see
#  <http://www.gnu.org/licenses/>.
#
######################################################################


# Remove insignificant volatile data from output files of tests
#
# Data affected:
#  timing data, date and time information in output files

# usage:
# normalize.pl <infile> <aspect_src_dir>

$infilename = shift;
$aspect_src_dir = shift;

open my $in,  '<', $infilename or die "Can't read file: $infilename!";

while (<$in>)
{
    # remove lines starting with --
    s/^--.*\n//;
    # remove lines starting and ending with | (from TimerOutput)
    s/^\|.*\|\n//;
    # replace lines starting with +-- and ending with --+ (from TimerOutput)
    s/^\+--.*--\+\n//;

    # remove time/date from gnuplot files:
    s/\# Time =.*//;
    s/\# Date =.*//;

    # replace current source directory by $ASPECT_SOURCE_DIR
    s/\Q$aspect_src_dir\E/\$ASPECT_SOURCE_DIR/g;

    # Exceptions
    s/line <\d+> of file <.*\//file </;

    # eat timestamp in graphical output
    s/by the deal.II library on.*//;

    print $_;
}
