set terminal pngcairo size 1600,1050 enhanced font 'Verdana,14'
set output 'figure_4.png'
set yrange [1e-14:1]
set xrange [0:50]
set format y '%g'
set logscale y
set multiplot title "Spiegelman benchmark convergence result for refinement 4 (64 by 16 cells)" font 'Verdana, 20'
set size 0.33,0.49
set origin 0,0.48
unset key
set title "Vel: 2.5cm/yr, \\eta_{ref}: 1e23, mLT: 9e-1"
plot \
'results/singleAdvectionanditeratedStokes_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'purple' t 'Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P150_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'purple' t 'DC Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P0_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'black' t '0 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P0_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'black' t '0 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P5_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'blue' t '5 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P5_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'blue' t '5 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P15_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'orange' t '15 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P15_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'orange' t '15 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P25_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'red' t '25 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P25_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'red' t '25 Picard stabilized'

set origin 0.33,0.48
set title "Vel: 5cm/yr, \\eta_{ref}: 1e24, mLT: 9e-1"
plot \
'results/singleAdvectionanditeratedStokes_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'purple' t 'Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P150_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'purple' t 'DC Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P0_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'black' t '0 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P0_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'black' t '0 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P5_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'blue' t '5 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P5_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'blue' t '5 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P15_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'orange' t '15 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P15_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'orange' t '15 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P25_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'red' t '25 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P25_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'red' t '25 Picard stabilized'

set origin 0.66,0.48
set title "Vel: 12.5cm/yr, \\eta_{ref}: 5e24, mLT: 9e-1"
plot \
'results/singleAdvectionanditeratedStokes_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'purple' t 'Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P150_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'purple' t 'DC Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P0_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'black' t '0 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P0_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'black' t '0 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P5_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'blue' t '5 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P5_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'blue' t '5 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P15_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'orange' t '15 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P15_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'orange' t '15 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_9e-1_P25_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'red' t '25 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_9e-1_P25_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'red' t '25 Picard stabilized'

set origin 0.0,0.0
set title "Vel: 2.5cm/yr, \\eta_{ref}: 1e23, mLT: 1e-8"
plot \
'results/singleAdvectionanditeratedStokes_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'purple' t 'Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P150_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'purple' t 'DC Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P0_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'black' t '0 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P0_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'black' t '0 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P5_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'blue' t '5 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P5_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'blue' t '5 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P15_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'orange' t '15 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P15_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'orange' t '15 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P25_vel_25_BV_1e23/plot.dat' u 9:10 w l lw 4 lc rgb 'red' t '25 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P25_vel_25_BV_1e23/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'red' t '25 Picard stabilized'

set origin 0.33,0.0
set title "Vel: 5cm/yr, \\eta_{ref}: 1e24, mLT: 1e-8"
plot \
'results/singleAdvectionanditeratedStokes_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'purple' t 'Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P150_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'purple' t 'DC Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P0_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'black' t '0 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P0_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'black' t '0 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P5_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'blue' t '5 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P5_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'blue' t '5 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P15_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'orange' t '15 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P15_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'orange' t '15 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P25_vel_50_BV_1e24/plot.dat' u 9:10 w l lw 4 lc rgb 'red' t '25 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P25_vel_50_BV_1e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'red' t '25 Picard stabilized'

set key bottom right
set key spacing 0.8
set origin 0.66,0.0
set title "Vel: 12.5cm/yr, \\eta_{ref}: 5e24, mLT: 1e-8"
plot \
'results/singleAdvectionanditeratedStokes_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'purple' t 'Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P150_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'purple' t 'DC Picard', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P0_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'black' t '0 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P0_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'black' t '0 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P5_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'blue' t '5 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P5_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'blue' t '5 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P15_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'orange' t '15 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P15_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'orange' t '15 Picard stabilized', \
'results/singleAdvectionandNewtonStokes_NS_none_mLT_1e-8_P25_vel_125_BV_5e24/plot.dat' u 9:10 w l lw 4 lc rgb 'red' t '25 Picard unstabilized', \
'results/singleAdvectionandNewtonStokes_NS_SPD_mLT_1e-8_P25_vel_125_BV_5e24/plot.dat' u 9:10 w p pt 1 lw 4 lc rgb 'red' t '25 Picard stabilized'
