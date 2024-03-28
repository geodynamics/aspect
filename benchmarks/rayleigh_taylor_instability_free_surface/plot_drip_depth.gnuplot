reset
set term png enhanced #dl 0.1 
set encoding utf8

set border lw 3
set size 0.7,1

set title "Drip depth over time for Rayleigh-Taylor instability with free surface"
set xlabel 'Time [kyr]' offset 0.5,0,0 font "Times-New-Roman:Bold, 12" 
set ylabel 'Drip depth [km]' offset 0,0.5,0 font "Times-New-Roman:Bold, 12" 
set output 'drip_depth.png'

set style line 1 lc rgb 'black' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 2 lc rgb 'brown' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 3 lc rgb '#D2691E' lt 1 lw 2 pt 7 #pi -1 ps 0.8 
set style line 4 lc rgb 'orange' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 5 lc rgb 'black' lt 1 lw 1 pt 7 pi -1 ps 0.8 
set style line 7 lc rgb 'yellow' lt 1 lw 2 pt 7 pi -1 ps 0.8 
set style line 6 lc rgb 'blue' lt 1 lw 2 pt 7 pi -1 ps 0.8 

set xrange [0:6e3]
set yrange [480:80]

set key bottom left

plot 'output_nostab_2500/statistics' u ($2/1e3):($19) ls 1 w l t 'no stab, dt=2500', \
     'output_nostab_5000/statistics' u ($2/1e3):($19) ls 2 w l t 'no stab, dt=5000', \
     'output_stab_5000/statistics' u ($2/1e3):($19) ls 3 w l  t 'stab, dt=5000', \
     'output_stab_5000_normal/statistics' u ($2/1e3):($19) ls 4 w l t 'stab, dt=5000, normal', \
     'Kaus_2010.csv' u ($1*1e3):(-$2) ls 6 w l t 'Kaus ea 2010', \
     0 ls 5 t ""
