reset
set term png
set output 'Topography_benchmark_1.png'

set xrange [0:1]
# all
set yrange [0:0.1]
set xlabel 'X [-]'
set ylabel 'Topography [m]'

set title 'Surface topography for benchmark 1'

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

plot \
     '1_sine_zero_flux/topography.00060' u 1:3 ls 1 with linespoints t 'ASPECT - t60', \
     '1_sine_zero_flux/topography.00120' u 1:3 ls 2 with linespoints t 'ASPECT - t120', \
     '1_sine_zero_flux/topography.00180' u 1:3 ls 3 with linespoints t 'ASPECT - t180', \
     '1_sine_zero_flux/topography.00060' u 1:4 ls 4 with lines t 'Analytical - t60', \
     '1_sine_zero_flux/topography.00120' u 1:4 ls 5 with lines notitle, \
     '1_sine_zero_flux/topography.00180' u 1:4 ls 6 with lines notitle, \
