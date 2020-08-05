#!/usr/bin/perl
######################################################################
#
# Copyright (C) 2017
#
# Remove insignificant volatile data from output files of tests
#
# Data affected:
#  JobID line (containing date)
#  line number of exceptions
#  start and final residual in iterations
#  small doubles
######################################################################

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

    # replace current source directory by ASPECT_DIR
    s/$aspect_src_dir/ASPECT_DIR/g;

    # Exceptions
    s/line <\d+> of file <.*\//file </;

    # eat timestamp in graphical output
    s/by the deal.II library on.*//;

    print $_;
}

# Several date and time strings

#s/%%Creation Date:.*//;
#s/\"created\".*//;
#s/^\s+Time =.*//;
#s/^\s+Date =.*//;
#s/Time tag:.*//g;



# Make small exponentials zero

#s/-?\d?\.\d+e-[123456789]\d+/0.00/g;

# See if we have a -0.0... (not followed by any other digit) and replace it
# by the same number without the negative sign
#s/-0\.(0+)(?!\d)/0.\1/g;

# Residual values

#s/value.*//;
#s/with residual.*//;


# remove deal.II debug output
#s/^DEAL.*::_.*\n//g;
