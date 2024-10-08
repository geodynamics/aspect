#!/usr/bin/env perl

$filename=$ARGV[0];
while(<STDIN>)
{
    if ($filename eq "screen-output")
    {
	s/   Solving Stokes system (GMG)... (\d+)\+0 iterations./   Solving Stokes system (GMG)... XYZ iterations./;
    }
    print $_;
}
