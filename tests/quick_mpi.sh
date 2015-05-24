#!/usr/bin/env perl

$filename=$ARGV[0];
while(<STDIN>)
{
    if ($filename eq "screen-output")
    {
	s/   Solving Stokes system... (\d+) iterations./   Solving Stokes system... XYZ iterations./;
    }
    print $_;
}
