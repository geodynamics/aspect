# Plot results from Crameri et. al. case 1 benchmark
# Includes curves for Underworld, SULEC, MILAMIN_VEP, as well
# as an approximate analytical solution


plot "UNDERWORLD_fs.dat" using 1:2, \
     "MILAMIN_VEP.dat" using 1:2, \
     "SULEC_fs.dat" using 1:2, \
     "output/statistics" using ($2/1000):($15/1000), \
     "output/statistics" using ($2/1000):(7. * exp(-$2/14825)) with lines

pause -1

