set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
set output 'convection.png'

unset key
set view map
set logscale xy
set logscale cb
set cbrange [0.999:1.001]

set xlabel "Viscosity in Pa s"
set ylabel "Temperature variation across the mantle in K"

splot 'onset-convection-data.csv' using 1:2:($4/$3) with pm3d
