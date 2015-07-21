# Plot results from Crameri et. al. case 1 benchmark
# Includes curves for Underworld, SULEC, MILAMIN_VEP, as well
# as an approximate analytical solution


plot "UNDERWORLD_fs.dat" using 1:2, \
     "MILAMIN_fs.dat" using ($1/1000):2, \
     "SULEC_fs.dat" using 1:2, \
     "output/statistics" using ($2/1.e6):16

pause -1

