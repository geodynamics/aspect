# -------------------------------------------------------------------------------------
# GNUPLOT script to plot gravity acceleration and potential along the radius for
# a hollow sphere filled with a constant density. Are plotted both the numerical
# and analytical solution calculated in the gravity_point_values postprocessor. 
# Are also plotted the difference and error between gravity acceleration and potential 
# with their respective theoretical values.  
# -------------------------------------------------------------------------------------

set term pdf enhanced font "Times,8pt"
set grid

set xlabel 'r (km)' font ',10'
set pointsize 0.4 

set style rect fc lt -1 fs solid 0.15 noborder
set obj rect from 3480, graph 0 to 6371, graph 1

set key inside 
set key font ", 10"

set ylabel 'g' font ',10'
set output 'profile_gravity_const_g.pdf'
plot [:15000][]\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):10 w p ps .5 t 'numerical',\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):11 w l lw .5 lt -1 t 'analytical',\

set ylabel 'U' font ',10'
set key inside left 
set output 'profile_gravity_const_U.pdf'
plot [:15000][]\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):12 w p ps .5 t 'numerical',\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):13 w l lw .5 lt -1 t 'analytical',\



set obj rect from 6596, graph 0 to 6621, graph 1
set key inside right 

set output 'profile_gravity_const_g_error.pdf'
plot [6371:15000][]\
 1e-5 lt -1 t '1 mGal',\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):(abs($10-$11)) w lp t 'numerical'\

set ylabel '|U_{meas}-U_{th}|' font ',10'
set output 'profile_gravity_const_U_error.pdf'
plot [6371:15000][]\
 1e-5 lt -1 t '1 mGal',\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):(abs($12-$13)) w lp t 'numerical'\


set log y
set format y "10^{%L}" 
set key outside

set ylabel '|g_{meas}-g_{th}|' font ',10'
set output 'profile_gravity_const_g_errorlog.pdf'
plot [6371:15000][]\
 1e-5 lt -1 t '1 mGal',\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):(abs($10-$11)) w lp t 'numerical'\

set ylabel '|U_{meas}-U_{th}|' font ',10'
set output 'profile_gravity_const_U_errorlog.pdf'
plot [6371:15000][]\
 1e-5 lt -1 t '1 mGal',\
 'output-gravity_constant_profile/output_gravity/gravity-00000' u ($1/1e3):(abs($12-$13)) w lp t 'numerical'\


