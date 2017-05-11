#!/usr/bin/env perl

$filename=$ARGV[0];
while(<STDIN>)
{
    if ($filename eq "screen-output")
    {
	s/   Solving Stokes system... 0\+(\d+) iterations./   Solving Stokes system... 0+XYZ iterations./;
    }
    print $_;
}
