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
4.87452472e+00 4.28499932e+01 case1a_ref3.stat
4.88459738e+00 4.28656773e+01 case1a_ref4.stat
4.88442981e+00 4.28650353e+01 case1a_ref5.stat
4.88441082e+00 4.28649626e+01 case1a_ref6.stat
4.88440900e+00 4.28649470e+01 case1a_reference.stat
8.97241383e+00 2.03825007e+02 case1b_ref3.stat
1.04949306e+01 1.93087841e+02 case1b_ref4.stat
1.05339127e+01 1.93214560e+02 case1b_ref5.stat
1.05339404e+01 1.93214792e+02 case1b_ref6.stat
1.05340950e+01 1.93214540e+02 case1b_reference.stat
1.65309062e+01 1.71646325e+03 case1c_ref3.stat
1.79259952e+01 8.71596864e+02 case1c_ref4.stat
2.18584327e+01 8.33321024e+02 case1c_ref5.stat
2.19710401e+01 8.33972141e+02 case1c_ref6.stat
2.19724650e+01 8.33989770e+02 case1c_reference.stat
5.49598094e+00 1.20550074e+03 case2a_ref3.stat
1.04029665e+01 5.46821574e+02 case2a_ref4.stat
1.01930353e+01 4.75819129e+02 case2a_ref5.stat
1.00711673e+01 4.80164721e+02 case2a_ref6.stat
1.00660000e+01 4.80433400e+02 case2a_reference.stat
6.34755300e+00 1.69487825e+02 case2b_ref3.stat
6.94577830e+00 1.70006734e+02 case2b_ref4.stat
6.92990946e+00 1.71555149e+02 case2b_ref5.stat
6.92964919e+00 1.71744538e+02 case2b_ref6.stat
6.92990000e+00 1.71755000e+02 case2b_reference.stat

steps needed:
72 output-case1a_ref3/statistics
128 output-case1a_ref4/statistics
232 output-case1a_ref5/statistics
470 output-case1a_ref6/statistics
699 output-case1b_ref3/statistics
344 output-case1b_ref4/statistics
690 output-case1b_ref5/statistics
1368 output-case1b_ref6/statistics
1231 output-case1c_ref3/statistics
1823 output-case1c_ref4/statistics
3399 output-case1c_ref5/statistics
6726 output-case1c_ref6/statistics
1482 output-case2a_ref3/statistics
3141 output-case2a_ref4/statistics
6641 output-case2a_ref5/statistics
13593 output-case2a_ref6/statistics
2720 output-case2b_ref3/statistics
5518 output-case2b_ref4/statistics
17034 output-case2b_ref5/statistics
34166 output-case2b_ref6/statistics
