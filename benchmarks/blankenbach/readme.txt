Blankenback benchmark

This folder allows running Blankenbach case 1a, 2b, 2a, and 2b from the paper

Blankenbach, B., et al. "A benchmark comparison for mantle convection codes."
Geophysical Journal International 98.1 (1989): 23-38.

The reference values (given in caseXX_reference.stat) are from that paper (see
Table 9). When you run the benchmarks by typing in "make" after building the
plugin in the plugin/ subfolder, you will be eventually presented with the
output below. Pleas note that these computations take a long time to reach
steady state, especially on finer meshes.

See the comment in the "Makefile" in this folder for more options for running
the computations.

Output:
# Nu           Vrms           name:
4.78661864e+00 4.34590432e+01 case1a_ref4.stat
4.87927972e+00 4.29377468e+01 case1a_ref5.stat
4.88993106e+00 4.28733838e+01 case1a_ref6.stat
4.88680525e+00 4.28659548e+01 case1a_ref7.stat
4.88440900e+00 4.28649470e+01 case1a_reference.stat
1.01134412e+01 2.09668282e+02 case1b_ref4.stat
1.02903266e+01 1.96513188e+02 case1b_ref5.stat
1.05150372e+01 1.93631922e+02 case1b_ref6.stat
1.05471266e+01 1.93263017e+02 case1b_ref7.stat
1.05340950e+01 1.93214540e+02 case1b_reference.stat
1.69902453e+01 1.01175072e+03 case1c_ref4.stat
2.11241960e+01 9.03094212e+02 case1c_ref5.stat
2.13365914e+01 8.49967378e+02 case1c_ref6.stat
2.18874917e+01 8.36119931e+02 case1c_ref7.stat
2.19724650e+01 8.33989770e+02 case1c_reference.stat
9.52888539e+00 6.83450568e+02 case2a_ref4.stat
9.96931688e+00 5.26034725e+02 case2a_ref5.stat
1.00761750e+01 4.91233242e+02 case2a_ref6.stat
1.00749710e+01 4.81948381e+02 case2a_ref7.stat
1.00660000e+01 4.80433400e+02 case2a_reference.stat
6.02809048e+00 1.78813766e+02 case2b_ref4.stat
6.91232842e+00 1.71856143e+02 case2b_ref5.stat
6.91232842e+00 1.71856143e+02 case2b_ref6.stat
6.92990000e+00 1.71755000e+02 case2b_reference.stat

steps needed:
     165 output-case1a_ref4/statistics
     281 output-case1a_ref5/statistics
     500 output-case1a_ref6/statistics
     968 output-case1a_ref7/statistics
     456 output-case1b_ref4/statistics
     794 output-case1b_ref5/statistics
    1452 output-case1b_ref6/statistics
    2778 output-case1b_ref7/statistics
    1971 output-case1c_ref4/statistics
    3684 output-case1c_ref5/statistics
    6899 output-case1c_ref6/statistics
   13501 output-case1c_ref7/statistics
    2627 output-case2a_ref4/statistics
    7335 output-case2a_ref5/statistics
   14711 output-case2a_ref6/statistics
   27397 output-case2a_ref7/statistics
    2770 output-case2b_ref4/statistics
    5579 output-case2b_ref5/statistics
   11190 output-case2b_ref6/statistics
    1908 output-case2b_ref7/statistics
