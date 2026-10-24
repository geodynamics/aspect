#!/usr/bin/perl
# ------------------------------------------------------
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# SPDX-FileCopyrightText: Copyright (C) 2017-2026 by the ASPECT authors.
#
# This file is part of ASPECT.
#
# Detailed license information governing the source code and
# contributions can be found in the folder LICENSES and in
# LICENSE and CONTRIBUTING.md at the top level directory.
#
# ------------------------------------------------------


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
