# -------------------------------------------------------------------------------------
# GNUPLOT script to plot gravity acceleration and potential at every query points for
# a sphere filled with the PREM density. We show the two gravity outputs from the two
# sampling scheme Fibonacci spiral and equiangular map, and their respective average 
# gravity acceleration or potential value.
# The average values are taken from the statsitics file in the output directory.
# -------------------------------------------------------------------------------------

set term pdf enhanced font "Times,8pt"
set grid

set xlabel 'index points' font ',10'
set pointsize 0.4 

set key outside 
set key font ", 10"

set ylabel 'g' font ',10'
set output 'stats_gravity_prem_g.pdf'
plot [0:2702][]\
 'output-gravity_prem_map/output_gravity/gravity-00000' u 10 w p ps .2 pt 1 lc rgb '#aadc32' t 'equiangular map',\
 'output-gravity_prem_spiral/output_gravity/gravity-00000' u 10 w p ps .2 pt 7 t 'Fibonacci spiral',\
 9.169056757829 w l lw 2 lc rgb '#aadc32' t 'mean g map',\
 9.168757879666 w l lw 2 lt 2 t 'mean g spiral'\

set ylabel 'U' font ',10'
set output 'stats_gravity_prem_U.pdf'
plot [0:2702][]\
 'output-gravity_prem_map/output_gravity/gravity-00000' u 12 w p ps .2 pt 1 lc rgb '#aadc32' t 'equiangular map',\
 'output-gravity_prem_spiral/output_gravity/gravity-00000' u 12 w p ps .2 pt 7 t 'Fibonacci spiral',\
 -6.047807249200e+07 w l lw 2 lc rgb '#aadc32' t 'mean U map',\
 -6.047713222550e+07 w l lw 2 lt 2 t 'mean U spiral'\

