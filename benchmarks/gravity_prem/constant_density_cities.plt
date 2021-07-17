# -------------------------------------------------------------------------------------
# GNUPLOT script to plot gravity acceleration and potential at 10 locations 225 km 
# above the surface of a hollow sphere filled with a constant density. We varied 
# the number of the quadrature degree increase from 2 to 4 to see a convergence of the 
# numerical solution towards the analytical solution also calculated in the
# gravity_point_values postprocessor.
# Are also plotted the difference and error between gravity acceleration 
# and potential with their respective theoretical values.  
# -------------------------------------------------------------------------------------

set term pdf enhanced font "Times,8pt"
set grid

set xlabel 'city index' font ',10'
set pointsize 0.4 

set style rect fc lt -1 fs solid 0.15 noborder
set obj rect from 3480, graph 0 to 6371, graph 1

set key inside 
set key font ", 10"

set ylabel 'g' font ',10'
set output 'cities_gravity_const_g.pdf'
plot [][]\
 'output-gravity_constant_list_nq2/output_gravity/gravity-00000' u 10 w lp ps .5 t 'numerical nq+2',\
 'output-gravity_constant_list_nq3/output_gravity/gravity-00000' u 10 w lp ps .5 t 'numerical nq+3',\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u 10 w lp ps .5 t 'numerical nq+4',\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u 11 w l lw .5 lt -1 t 'analytical',\

set ylabel 'U' font ',10'
set output 'cities_gravity_const_U.pdf'
plot [][]\
 'output-gravity_constant_list_nq2/output_gravity/gravity-00000' u 12 w lp ps .5 t 'numerical nq+2',\
 'output-gravity_constant_list_nq3/output_gravity/gravity-00000' u 12 w lp ps .5 t 'numerical nq+3',\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u 12 w lp ps .5 t 'numerical nq+4',\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u 13 w l lw .5 lt -1 t 'analytical',\


unset key
set obj rect from 6596, graph 0 to 6621, graph 1

set ylabel '|g_{meas}-g_{th})|' font ',10'
set output 'cities_gravity_constant_g_error.pdf'
plot [][]\
 'output-gravity_constant_list_nq2/output_gravity/gravity-00000' u (abs($10-$11)) w lp ps 1.5,\
 'output-gravity_constant_list_nq3/output_gravity/gravity-00000' u (abs($10-$11)) w lp ps 1.5,\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u (abs($10-$11)) w lp ps 1.5,\

set ylabel '|U_{meas}-U_{th}|' font ',10'
set output 'cities_gravity_const_U_error.pdf'
plot [][]\
 'output-gravity_constant_list_nq2/output_gravity/gravity-00000' u (abs($12-$13)) w lp ps 1.5,\
 'output-gravity_constant_list_nq3/output_gravity/gravity-00000' u (abs($12-$13)) w lp ps 1.5,\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u (abs($12-$13)) w lp ps 1.5,\

set log y
set format y "10^{%L}" 

set ylabel '|(g_{meas}-g_{th})/g_{th}|' font ',10'
set output 'cities_gravity_constant_g_errorlog.pdf'
plot [][]\
 'output-gravity_constant_list_nq2/output_gravity/gravity-00000' u (abs(($10-$11)/$11)) w lp ps 1.5,\
 'output-gravity_constant_list_nq3/output_gravity/gravity-00000' u (abs(($10-$11)/$11)) w lp ps 1.5,\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u (abs(($10-$11)/$11)) w lp ps 1.5,\

set log y
set ylabel '|(U_{meas}-U_{th})/U_{th}|' font ',10'
set output 'cities_gravity_constant_U_errorlog.pdf'
plot [][]\
 'output-gravity_constant_list_nq2/output_gravity/gravity-00000' u (abs(($12-$13)/$13)) w lp ps 1.5,\
 'output-gravity_constant_list_nq3/output_gravity/gravity-00000' u (abs(($12-$13)/$13)) w lp ps 1.5,\
 'output-gravity_constant_list_nq4/output_gravity/gravity-00000' u (abs(($12-$13)/$13)) w lp ps 1.5,\


