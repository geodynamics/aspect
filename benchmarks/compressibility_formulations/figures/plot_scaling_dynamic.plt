# Gnuplot script

set terminal pdf color solid dashed font "Arial,12" size 7.5cm, 24cm
set output 'mass_flux_error_dynamic.pdf'

set multiplot layout 3,1
set size ratio 1.0


set title "Mass flux error" font "Arial,14" 
set xlabel "#Timesteps" font "Arial,14" 
set ylabel "Mass flux error" font "Arial,14" 
set logscale xy
set format y "10^{%S}"
#set format x "10^{%S}"
set ytics (1,1e-3,1e-6,1e-9)


set xrange [0.8:160]
set yrange [1e-9:10]

set grid ytics
set key at 140,0.3 noautotitles

f(x)= 3e-4 / (x*x)
g(x)= 3e-3 / x

plot "../lateral-pipe-transient/mass_flux_error" using (31557600e7/$1):2 with linespoints linetype 1 pointtype 7 ps 1.2 linecolor rgb "#4b03a1" lw 3 title "Anelastic liquid approximation (ALA)", \
     "../lateral-pipe-transient/mass_flux_error" using (31557600e7/$1):3 with linespoints linetype 1 pointtype 9 ps 1.2 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \
     "../lateral-pipe-transient/mass_flux_error" using (31557600e7/$1):4 with linespoints linetype 1 pointtype 11 ps 1.2 linecolor rgb "#a82296" lw 3 title "Hydrostatic compression", \
     "../lateral-pipe-transient/mass_flux_error" using (31557600e7/$1):5 with linespoints linetype 1 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density", \
     f(x) title '2nd order' with lines linestyle 2 linecolor "gray" lw 3
     # "./mass_flux_error" using (31557600e7/$1):3 with linespoints linetype 1 pointtype 9 ps 0.7 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \

set title "Mass flux error" font "Arial,14" 
set xlabel "#Timesteps" font "Arial,14" 
set ylabel "Mass flux error" font "Arial,14" 
set format y "10^{%S}"
#set format x "10^{%S}"
set ytics (1,1e-3,1e-6,1e-9)
unset yrange
set xrange [0.8:160]
set yrange [1e-9:10]

set grid ytics
set key at 100,1e-5 noautotitles

f(x)= 1e-1 / (x*x)
g(x)= 3e-3 / x

plot "../lateral-pipe-increase-pressure/mass_flux_error" using (31557600e7/$1):2 with linespoints linetype 1 pointtype 7 ps 1.2 linecolor rgb "#4b03a1" lw 3 title "Anelastic liquid approximation", \
     "../lateral-pipe-increase-pressure/mass_flux_error" using (31557600e7/$1):3 with linespoints linetype 1 pointtype 9 ps 1.2 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \
     "../lateral-pipe-increase-pressure/mass_flux_error" using (31557600e7/$1):4 with linespoints linetype 1 pointtype 11 ps 1.2 linecolor rgb "#a82296" lw 3 title "Hydrostatic compression", \
     "../lateral-pipe-increase-pressure/mass_flux_error" using (31557600e7/$1):5 with linespoints linetype 1 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density", \
     "../lateral-pipe-increase-pressure/mass_flux_error" using (31557600e7/$1):6 with linespoints linetype 1 dt 2 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density (full pressure)", \
     f(x) title '2nd order' with lines linestyle 2 linecolor "gray" lw 3
     # "./mass_flux_error" using (31557600e7/$1):3 with linespoints linetype 1 pointtype 9 ps 0.7 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \



set title "Mass flux error" font "Arial,14" 
set xlabel "#Cells" font "Arial,14" 
set ylabel "Mass flux error" font "Arial,14" 
set logscale xy
set format y "10^{%S}"
#set format x "10^{%S}"
set ytics (1,1e-3,1e-6,1e-9,1e-12)


set xrange [0.8:160]
set yrange [1e-13:250]

set grid ytics
set key at 140,150 noautotitles

f(x)= 3e-4 / (x*x)
g(x)= 3e-3 / x

plot "../lateral-pipe-advect-smooth/mass_flux_error" using 1:2 with linespoints linetype 1 pointtype 7 ps 1.2 linecolor rgb "#4b03a1" lw 3 title "Anelastic liquid approximation", \
     "../lateral-pipe-advect-smooth/mass_flux_error" using 1:3 with linespoints linetype 1 pointtype 9 ps 1.2 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \
     "../lateral-pipe-advect-smooth/mass_flux_error" using 1:4 with linespoints linetype 1 pointtype 11 ps 1.2 linecolor rgb "#a82296" lw 3 title "Hydrostatic compression", \
     "../lateral-pipe-advect-smooth/mass_flux_error" using 1:5 with linespoints linetype 1 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density", \
     f(x) title '2nd order' with lines linestyle 2 linecolor "gray" lw 3
     # "./mass_flux_error" using (31557600e7/$1):3 with linespoints linetype 1 pointtype 9 ps 0.7 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \

unset multiplot
#replot

