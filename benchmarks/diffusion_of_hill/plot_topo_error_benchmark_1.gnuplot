reset
set term png
set output 'Topography_error_benchmark_1.png'

set xrange [0:1]
set yrange [0:2]
set xlabel 'X [-]'
set ylabel 'Error topography [%]'

set title 'Surface topography error for benchmark 1'

set style line 1 lw 2 lc rgb "black"
set style line 2 lw 3 lc rgb "black"
set style line 3 lw 4 lc rgb "black"

set style line 4 lw 2 lc rgb "dark-grey"
set style line 5 lw 3 lc rgb "dark-grey"
set style line 6 lw 4 lc rgb "dark-grey"

set style line 7 lw 2 lc rgb "light-grey"
set style line 8 lw 3 lc rgb "light-grey"
set style line 9 lw 4 lc rgb "light-grey"

set style line 10 lw 2 lc rgb "light-coral"
set style line 11 lw 3 lc rgb "light-coral"
set style line 12 lw 4 lc rgb "light-coral"

set key top center

round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
round2(x, n) = round(x*10**n)*10.0**(-n)

plot \
     '1_sine_zero_flux/topography.00060' u ($1):(abs((round2(($4),5)-$3)/round2(($4),5)*100.)) ls 1 with linespoints t 't60', \
     '1_sine_zero_flux/topography.00120' u ($1):(abs((round2(($4),5)-$3)/round2(($4),5)*100.)) ls 2 with linespoints t 't120', \
     '1_sine_zero_flux/topography.00180' u ($1):(abs((round2(($4),5)-$3)/round2(($4),5)*100.)) ls 3 with linespoints t 't180', \
