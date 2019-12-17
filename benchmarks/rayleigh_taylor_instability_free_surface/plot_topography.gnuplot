reset
set term png enhanced dl 0.1 
set encoding utf8

set border lw 3
set size 0.7,1

set title "Rayleigh-Taylor instability with free surface"
set xlabel 'Time [My]' offset 0.5,0,0 font "Times-New-Roman:Bold, 12" 
set ylabel 'Maximum topography [m]' offset 0,0.5,0 font "Times-New-Roman:Bold, 12" 
set output 'topography.png'

set style line 1 lc rgb 'black' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 2 lc rgb 'brown' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 3 lc rgb '#D2691E' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 4 lc rgb 'orange' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 5 lc rgb 'black' lt 1 lw 2 pt 7 pi -1 ps 0.8 

set xrange [0:6]
set yrange [0:2100]

set key bottom left

plot 'output_nostab_2500/statistics' u ($2/1e6):15 ls 1 t 'no stab., dt=2500 yr', \
     'output_nostab_10000/statistics' u ($2/1e6):15 ls 2 w l t 'no stab., dt=5000 yr', \
     'output_stab_10000/statistics' u ($2/1e6):15 ls 3 w l t 'stab., dt=5000 yr', \
     'output_stab_10000_normal/statistics' u ($2/1e6):15 ls 4 w l t 'stab., dt=5000 yr, normal', \
     'output_Rose_stab_CFL02/statistics' u ($2/1e6):15 ls 5 w l t 'stab., CFL 0.2, 1 km', \

