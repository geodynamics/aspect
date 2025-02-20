# The Rayleigh-Taylor instability benchmark

The shell script `run_benchmark` loops over viscosity values and wavelength values and
for each combination generates an input file `inputfile.prm` based on the `blank.prm` file
to which additional commands are added.
The code is then run with inputfile.prm and the standard output is redirected in a
temporary file 'opla' of no importance. Finally the viscosity and wavelength
values are written to the 'vy.dat' file followed by the last line of the statistics.
This file can later be processed with the 'gnuplot_script'.

There is also the standalone input file `rayleigh_taylor_instability.prm` which
runs only one of the cases (isoviscous and lambda=256km).

See [](doc/rayleigh_taylor_instability) for more information.
