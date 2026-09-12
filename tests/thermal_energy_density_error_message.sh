#!/usr/bin/env perl

# This script processes the output of the thermal energy density plugin failing test.
# Special processing is needed because the error information contains the whole path,
# which is not repoducable in the tester and the ASPECT tester outputs the dealii part 
# of the error message differently that what is ouput in the ASPECT configuration on 
# which the test was made. This is based on the default script in tests/cmake/default.

$filename=$ARGV[0];
my $to_find="An error occurred in";
my $to_replace="(line in output replaced by thermal_energy_density.sh script)";
my $skip_rest=0;
while(<STDIN>)
{
    if ($filename eq "screen-output")
    {
        next if $_ =~ m/^-- /;
        next if $_ =~ m/^\|/;
	    next if $skip_rest>0;

	    s/$to_find.+/$to_replace/;
	    if($_ =~ m/^Aborting/)
	    {
		    $skip_rest=1;
	    }
    }
    print $_;
}
