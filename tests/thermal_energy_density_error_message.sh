#!/usr/bin/env perl

$filename=$ARGV[0];
my $to_find1="An error occurred in line <";
my $to_find2="> of file ";
my $to_find3="/source/numerics/derivative_approximation.cc> in function";
my $to_replace1="An error occurred in line <";
my $to_replace2="> of file <DEALII_DIR/source/numerics/derivative_approximation.cc> in function\n(line in output replaced by thermal_energy_density.sh script)";
while(<STDIN>)
{
    if ($filename eq "screen-output")
    {
	    s/$to_find1(.+)$to_find2.+$to_find3/$to_replace1$1$to_replace2/;
    }
    print $_;
}
