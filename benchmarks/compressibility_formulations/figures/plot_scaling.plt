# Gnuplot script

set terminal pdf color solid dashed font "Arial,12" size 15cm, 18cm
set output 'mass_flux_error.pdf'

set multiplot layout 2,2
set size ratio 1.0


set title "Mass flux error: Adiabatic temperature" font "Arial,14" 
set xlabel "#Cells" font "Arial,14" 
set ylabel "Mass flux error" font "Arial,14" 
set logscale xy
set format y "10^{%S}"
set ytics (1,1e-3,1e-6,1e-9)

set xrange [0.8:160]
set yrange [1e-10:10]

set grid ytics
set key top right noautotitles

f(x)= 2e-4 / (x*x)

plot "../vertical_pipe/mass_flux_error_adiabatic" using 1:2 with linespoints linetype 1 pointtype 7 ps 1.2 linecolor rgb "#4b03a1" lw 3 title "Anelastic liquid approximation (ALA)", \
     "../vertical_pipe/mass_flux_error_adiabatic" using 1:3 with linespoints linetype 1 pointtype 9 ps 1.2 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \
     "../vertical_pipe/mass_flux_error_adiabatic" using 1:4 with linespoints linetype 1 pointtype 11 ps 1.2 linecolor rgb "#a82296" lw 3 title "Hydrostatic compression", \
     "../vertical_pipe/mass_flux_error_adiabatic" using 1:5 with linespoints linetype 1 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density", \
     f(x) title '2nd order' with lines linestyle 2 linecolor "gray" lw 3 



set size ratio 1.0
set key off

set title "Mass flux error: Sub-adiabatic temperature" font "Arial,14" 
set xlabel "#Cells" font "Arial,14" 
set ylabel "Mass flux error" font "Arial,14" 
set logscale xy

set xrange [0.8:160]
set yrange [1e-10:10]

set grid ytics

set key top right noautotitles

plot "../vertical_pipe/mass_flux_error_sub_adiabatic" using 1:2 with linespoints linetype 4 pointtype 7 ps 1.2 linecolor rgb "#4b03a1" lw 3 title "Anelastic liquid approximation (ALA)", \
     "../vertical_pipe/mass_flux_error_sub_adiabatic" using 1:3 with linespoints linetype 1 pointtype 9 ps 1.2 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \
     "../vertical_pipe/mass_flux_error_sub_adiabatic" using 1:4 with linespoints linetype 4 pointtype 11 ps 1.2 linecolor rgb "#a82296" lw 3 title "Hydrostatic compression", \
     "../vertical_pipe/mass_flux_error_sub_adiabatic" using 1:5 with linespoints linetype 4 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density", \
     f(x) title '2nd order' with lines linestyle 2 linecolor "gray" lw 3 

set title "Mass flux error: Adiabatic temperature" font "Arial,14" 
set xlabel "#Cells" font "Arial,14" 
set ylabel "Mass flux error" font "Arial,14" 
set logscale xy

set xrange [0.8:160]
set yrange [1e-9:10]

set grid ytics
set key top right noautotitles

f(x)= 2e-4 / (x*x)

plot "../lateral_pipe/mass_flux_error_adiabatic" using 1:2 with linespoints linetype 1 pointtype 7 ps 1.2 linecolor rgb "#4b03a1" lw 3 title "Anelastic liquid approximation (ALA)", \
     "../lateral_pipe/mass_flux_error_adiabatic" using 1:3 with linespoints linetype 1 pointtype 9 ps 1.2 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \
     "../lateral_pipe/mass_flux_error_adiabatic" using 1:4 with linespoints linetype 1 pointtype 11 ps 1.2 linecolor rgb "#a82296" lw 3 title "Hydrostatic compression", \
     "../lateral_pipe/mass_flux_error_adiabatic" using 1:5 with linespoints linetype 1 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density", \
     f(x) title '2nd order' with lines linestyle 2 linecolor "gray" lw 3 



set size ratio 1.0
set key off

set title "Mass flux error: Sub-adiabatic temperature" font "Arial,14" 
set xlabel "#Cells" font "Arial,14" 
set ylabel "Mass flux error" font "Arial,14" 
set logscale xy

set xrange [0.8:160]
set yrange [1e-9:10]

set grid ytics

set key top right noautotitles

plot "../lateral_pipe/mass_flux_error_sub_adiabatic" using 1:2 with linespoints linetype 4 pointtype 7 ps 1.2 linecolor rgb "#4b03a1" lw 3 title "Anelastic liquid approximation (ALA)", \
     "../lateral_pipe/mass_flux_error_sub_adiabatic" using 1:3 with linespoints linetype 1 pointtype 9 ps 1.2 linecolor rgb "#fdc328" lw 3 title "Isentropic compression", \
     "../lateral_pipe/mass_flux_error_sub_adiabatic" using 1:4 with linespoints linetype 4 pointtype 11 ps 1.2 linecolor rgb "#a82296" lw 3 title "Hydrostatic compression", \
     "../lateral_pipe/mass_flux_error_sub_adiabatic" using 1:5 with linespoints linetype 4 pointtype 13 ps 1.2 linecolor rgb "#e56b5d" lw 3 title "Projected density", \
     f(x) title '2nd order' with lines linestyle 2 linecolor "gray" lw 3 

unset multiplot
