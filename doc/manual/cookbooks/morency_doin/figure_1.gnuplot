set terminal pdf size 6,3 enhanced color
set output "morency_doin_2004_fig1.pdf"

set ylabel "Depth (km)"
set yrange [-750:0]
unset key
set multiplot
set grid

set size 0.33,1.0
set origin 0.0,0.0
set title "Temperature"
set xlabel "Temperature (K)"
set xtics 0,1000,2000
set xrange [0:2000]
plot "< grep -v -e \"^$\" -e \" 0 $\" output/depth_average.gnuplot" using ($3-293):(-$1/1000) with lines

unset ylabel
set size 0.66,1.0
set origin 0.33,0.0
set title "Viscosity"
set xlabel "Viscosity (Pa s)"
set logscale x
set xrange [1e18:1e28]
set xtics 1e18,100,1e30
plot "< grep -v -e \"^$\" -e \" 0 $\" output/depth_average.gnuplot" using 11:(-$1/1000) with lines
unset logscale x
