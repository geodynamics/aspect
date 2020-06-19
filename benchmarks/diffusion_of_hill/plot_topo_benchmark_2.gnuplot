reset
set term png
set output 'Topography_benchmark_2.png'

set xrange [0:10e3]
# all
set yrange [0:110]
set xlabel 'X [-]'
set ylabel 'Topography [m]'

set title 'Surface topography for benchmark 2'

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

set key bottom center

plot \
     '2_sine_constant_h/topography.00030' u 1:3 ls 1 with linespoints t 'ASPECT - t30', \
     '2_sine_constant_h/topography.00060' u 1:3 ls 2 with linespoints t 'ASPECT - t60', \
     '2_sine_constant_h/topography.00090' u 1:3 ls 3 with linespoints t 'ASPECT - t90', \
     '2_sine_constant_h/topography.00030' u 1:4 ls 4 with lines t 'Analytical - t30', \
     '2_sine_constant_h/topography.00060' u 1:4 ls 5 with lines notitle, \
     '2_sine_constant_h/topography.00090' u 1:4 ls 6 with lines notitle, \
