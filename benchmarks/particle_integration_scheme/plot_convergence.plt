set terminal svg size 1500,900 font "Arial,20" enhanced background rgb 'white'
set output 'convergence.svg'

set multiplot layout 1,2

# Plot A
set size ratio 1.0
set origin 0.0,0.0
set lmargin at screen 0.1
set rmargin at screen 0.45
set tmargin at screen 0.9
set bmargin at screen 0.25

set xrange[2.5e-2:1.2]
set yrange[1e-14:1.0]
set xlabel "CFL Number" font "Arial,18"
set ylabel "L2 errors" font "Arial,18"
set logscale xy

f(x) = 0.1 * x
g(x) = 0.005 * x*x
h(x) = 1e-6 * x*x*x*x
set key font "Arial,16" outside center bottom samplen 1 maxrows 3 box

set label 1 "A)" at screen 0.25, screen 0.95 font "Arial,22"

plot "spatial/max_deviation_euler" using 1:2 with linespoints linetype 1 lw 3 pointtype 7 ps 1.0 linecolor rgb "#0A5F02" title "Forward Euler", \
     "spatial/max_deviation_rk2" using 1:2 with linespoints linetype 1 lw 3 pointtype 9 ps 1.0 linecolor rgb "#1c567a" title "RK2", \
     "spatial/max_deviation_rk4" using 1:2 with linespoints linetype 1 lw 3 pointtype 11 ps 1.0 linecolor rgb "#814292" title "RK4", \
     f(x) title 'linear' with lines dashtype 2 lw 3 linecolor rgb "gray", \
     g(x) title 'quadratic' with lines dashtype 3 lw 3 linecolor rgb "gray", \
     h(x) title 'quartic' with lines dashtype 4 lw 3 linecolor rgb "gray"

# Plot B
set size ratio 1.0
set origin 0.5,0.0
set lmargin at screen 0.55
set rmargin at screen 0.9
set tmargin at screen 0.9
set bmargin at screen 0.25

set xrange[2.5e-2:1.2]
set yrange[1e-14:1.0]
set xlabel "CFL Number" font "Arial,18"
set ylabel "L2 errors" font "Arial,18"
set logscale xy

f(x) = 0.05 * x
g(x) = 2e-3 * x*x
h(x) = 1e-7 * x*x*x*x
set key font "Arial,16" outside center bottom samplen 1 maxrows 4 box

set label 2 "B)" at screen 0.75, screen 0.95 font "Arial,22"

plot "time_exponential/max_deviation_euler" using 1:2 with linespoints linetype 1 lw 3 pointtype 7 ps 1.0 linecolor rgb "#0A5F02" notitle , \
     "time_exponential/max_deviation_rk2" using 1:2 with linespoints linetype 1 lw 3 pointtype 8 ps 1.0 linecolor rgb "#1c567a" notitle , \
     "time_exponential/max_deviation_rk4" using 1:2 with linespoints linetype 1 lw 3 pointtype 11 ps 1.0 linecolor rgb "#814292" notitle, \
     "time_exponential/max_deviation_rk4_analytic" using 1:2 with linespoints linetype 2 lw 3 pointtype 13 ps 1.0 linecolor rgb "#814292" title "RK4 with analytic halfstep", \
     f(x) notitle with lines dashtype 2 lw 3 linecolor rgb "gray", \
     g(x) notitle with lines dashtype 3 lw 3 linecolor rgb "gray", \
     h(x) notitle with lines dashtype 4 lw 3 linecolor rgb "gray"

unset multiplot
