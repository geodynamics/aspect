#!/usr/bin/env perl

$filename=$ARGV[0];
my $to_find="An error occurred in";
my $to_replace="(line in output replaced by thermal_energy_density.sh script)";
my $skip_rest=0;
while(<STDIN>)
{
    if ($filename eq "screen-output")
    {
<<<<<<< HEAD
	    s/$to_find1(.+)$to_find2.+$to_find3/$to_replace1$1$to_replace2/;
=======
            next if $_ =~ m/^-- /;
            next if $_ =~ m/^\|/;
	    next if $skip_rest>0;

	    s/$to_find.+/$to_replace/;
	    if($_ =~ m/^Aborting.*/)
	    {
		$skip_rest=1;
	    }
>>>>>>> 60126fc4b... fix script and test result for difference in dealii output.
    }
    print $_;
}
