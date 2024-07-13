set terminal qt size 1500,1000 enhanced font 'Verdana,14'
set format y '%g'
set yrange [1e-6:1]
set xrange [0:50]
set log y
set multiplot title "Tosi case 4 benchmark convergence results" font 'Verdana, 20' #No stabilisation of velocity block with a stress exponent of 3 and pressure boundary conditions" font 'Verdana, 20'
set size 0.5,0.49
set origin 0,0.48
unset key
set title "0-50 iterations, maximum LT: 9e-1"
plot \
'results/a9TS_itAdandSt_ST1e-20_UFSfalse_NSP-none_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-8_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:10 w l lw 2 t 'Picard', \
'results/a9TS_NS_ST1e-20_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Defect correction Picard', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton without stabilsation', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-SPD_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton with stabilsation'
set origin 0,0
set xrange [500:600]
set title "500-600 iterations, maximum LT: 9e-1"
plot \
'results/a9TS_itAdandSt_ST1e-20_UFSfalse_NSP-none_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-8_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:10 w l lw 2 t 'Picard', \
'results/a9TS_NS_ST1e-20_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Defect correction Picard', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton without stabilsation', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-SPD_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton with stabilsation'
set key top font 'Verdana,13'
set origin 0.5,0.48
set xrange [0:50]
set title "0-50 iterations, maximum LT: 1e-2"
plot \
'results/a9TS_itAdandSt_ST1e-20_UFSfalse_NSP-none_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-8_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:10 w l lw 2 t 'Picard', \
'results/a9TS_NS_ST1e-20_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Defect correction Picard', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT1e-2_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton without stabilsation', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-SPD_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT1e-2_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton with stabilsation'
unset key
set origin 0.5,0
set xrange [500:600]
set title "500-600 iterations, maximum LT: 1e-2"
plot \
'results/a9TS_itAdandSt_ST1e-20_UFSfalse_NSP-none_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-8_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:10 w l lw 2 t 'Picard', \
'results/a9TS_NS_ST1e-20_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT9e-1_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Defect correction Picard', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-none_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT1e-2_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton without stabilsation', \
'results/a9TS_NS_ST1e-2_UFSfalse_NSP-SPD_NSA-SPD_C0_g6_ag0_AEWfalse_UDSfalse_SF9.99e-1_NLT1e-5_ABT1e-2_LT1e-5_mLT1e-2_I150_P150_EW1_theta1_LS0_RSMfalse_AV-1/plot.dat' u 0:11 w l lw 2 t 'Newton with stabilsation'
pause mouse
