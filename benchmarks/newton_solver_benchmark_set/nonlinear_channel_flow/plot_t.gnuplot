set terminal png size 1500,1000 enhanced font 'Verdana,14'
set output 'figure_t.png'
set yrange [1e-14:1]
set xrange [0:60]
set format y '%g'
set logscale y
set multiplot title "With a minimum linear tolerance of 1e-8, a stress exponent of 3 and pressure boundary conditions" font 'Verdana, 20'
set size 0.51,0.49
set origin 0,0.48
unset key
set title "Without stabilization of the velocity block and without the RSM"
plot \
'results/BT_t_singleAdvectioniteratedStokes/plot.dat' u 9:10 w l lw 2 lc rgb 'purple' t 'Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_150_LS_0_RSM_false/plot.dat' u 9:10 w l lw 2 lc rgb 'green' t 'DC Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_0_RSM_false/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'blue' t 'Newton solver, 5 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 7 t 'Newton solver, 5 Picard, line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_0_RSM_false/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'orange' t 'Newton solver, 10 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 7 ps 1 t 'Newton solver, 10 Picard, line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_0_RSM_false/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'red' t 'Newton solver, 15 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 7 ps 1 t 'Newton solver, 15 Picard, line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_150_LS_0_RSM_false/plot.dat' u 9:10 w p lw 2 lc rgb 'green' pt 5 t 'DC Picard', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_0_RSM_false/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'blue' pt 5 t 'Newton solver, 5 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 5 t 'Newton solver, 5 Picard, line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_0_RSM_false/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'orange' pt 5 t 'Newton solver, 10 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 5 ps 1 t 'Newton solver, 10 Picard, with line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_0_RSM_false/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'red' pt 5 t 'Newton solver, 15 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 5 ps 1 t 'Newton solver, 15 Picard, with line search'
set origin 0.49,0.48

set key spacing 0.8
set key
set title "With stabilization of the velocity block and without the RSM"
plot \
'results/BT_t_singleAdvectioniteratedStokes/plot.dat' u 9:10 w l lw 2 lc rgb 'purple' t 'Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_150_LS_0_RSM_false/plot.dat' u 9:10 w l lw 2 lc rgb 'green' t 'DC Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_0_RSM_false/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'blue' t 'Newton solver, 5 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 7 t 'Newton solver, 5 Picard, line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_0_RSM_false/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'orange' t 'Newton solver, 10 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 7 ps 1 t 'Newton solver, 10 Picard, with line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_0_RSM_false/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'red' t 'Newton solver, 15 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 7 ps 1 t 'Newton solver, 15 Picard, with line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_150_LS_0_RSM_false/plot.dat' u 9:10 w p lw 2 lc rgb 'green' pt 5 t '', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_0_RSM_false/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'blue' pt 5 t '', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 5 t '', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_0_RSM_false/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'orange' pt 5 t '', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 5 ps 1 t '', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_0_RSM_false/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'red' pt 5 t '', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_100_RSM_false/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 5 ps 1 t ''
unset key
set origin 0,0
set title "Without stabilization of the velocity block and with the RSM"
plot \
'results/BT_t_singleAdvectioniteratedStokes/plot.dat' u 9:10 w l lw 2 lc rgb 'purple' t 'Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_150_LS_0_RSM_true/plot.dat' u 9:10 w l lw 2 lc rgb 'green' t 'DC Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_0_RSM_true/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'blue' t 'Newton solver, 5 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 7 t 'Newton solver, 5 Picard, line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_0_RSM_true/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'orange' t 'Newton solver, 10 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 7 ps 1 t 'Newton solver, 10 Picard, with line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_0_RSM_true/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'red' t 'Newton solver, 15 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 7 ps 1 t 'Newton solver, 15 Picard, with line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_150_LS_0_RSM_true/plot.dat' u 9:10 w p lw 2 lc rgb 'green' pt 5 t 'DC Picard', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_0_RSM_true/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'blue' pt 5 t 'Newton solver, 5 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_5_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 5 t 'Newton solver, 5 Picard, line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_0_RSM_true/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'orange' pt 5 t 'Newton solver, 10 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_10_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 5 ps 1 t 'Newton solver, 10 Picard, with line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_0_RSM_true/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'red' pt 5 t 'Newton solver, 15 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_none_mLT_1e-8_P_15_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 5 ps 1 t 'Newton solver, 15 Picard, with line search'
set origin 0.49,0
set title "With stabilization of the velocity block and with the RSM"
plot \
'results/BT_t_singleAdvectioniteratedStokes/plot.dat' u 9:10 w l lw 2 lc rgb 'purple' t 'Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_150_LS_0_RSM_true/plot.dat' u 9:10 w l lw 2 lc rgb 'green' t 'DC Picard', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_0_RSM_true/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'blue' t 'Newton solver, 5 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 7 t 'Newton solver, 5 Picard, line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_0_RSM_true/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'orange' t 'Newton solver, 10 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 7 ps 1 t 'Newton solver, 10 Picard, with line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_0_RSM_true/plot.dat' u 9:10 w l dt 1 lw 2 lc rgb 'red' t 'Newton solver, 15 Picard, no line search', \
'results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 7 ps 1 t 'Newton solver, 15 Picard, with line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_150_LS_0_RSM_true/plot.dat' u 9:10 w p lw 2 lc rgb 'green' pt 5 t 'DC Picard', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_0_RSM_true/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'blue' pt 5 t 'Newton solver, 5 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_5_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'blue' pt 5 t 'Newton solver, 5 Picard, line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_0_RSM_true/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'orange' pt 5 t 'Newton solver, 10 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_10_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'orange' pt 5 ps 1 t 'Newton solver, 10 Picard, with line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_0_RSM_true/plot.dat' u 9:10 w p dt 1 lw 2 lc rgb 'red' pt 5 t 'Newton solver, 15 Picard, no line search', \
'<tail -n 1 results/BT_t_singleAdvectioniteratedNewtonStokes_stabilization_SPD_mLT_1e-8_P_15_LS_100_RSM_true/plot.dat' u 9:10 w p dt 3 lw 2 lc rgb 'red' pt 5 ps 1 t 'Newton solver, 15 Picard, with line search'
