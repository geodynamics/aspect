#!/bin/gnuplot -persist

set term qt 0 font "Sans, 16"
set log x
set log y
set xrange [ 5e3 : 3e8 ]
set xlabel "Time [yr]"
set yrange [ 1e-7 : 1e-2 ]
set ylabel "Grain size [m]"
set format x "%2.0tx10^{%L}"
set format y "%2.0tx10^{%L}"
set grid mxtics mytics
plot "output-grain_size_plunge/statistics" using 2:14 with lines linecolor "green" linewidth 2 title "ASPECT Grain size [m]"

set term qt 1 font "Sans, 16"
set log x
set log y
set xrange [ 5e3 : 3e8 ]
set xlabel "Time [yr]"
set yrange [ 1e18 : 1e23 ]
set ylabel "Viscosity [Pa s]"
set grid mxtics mytics
plot "output-grain_size_plunge/statistics" using 2:25 with lines linecolor "green" linewidth 2 title "ASPECT Viscosity [m]"