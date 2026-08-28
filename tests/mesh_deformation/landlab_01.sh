#!/usr/bin/env perl

$filename=$ARGV[0];
while(<STDIN>)
{
    if ($filename eq "screen-output")
    {
        next if /Matplotlib is building the font cache/;
    }
    print $_;
}
