reset
set term png enhanced
set encoding utf8

set border lw 3
set size 0.7,1

set title "Rayleigh-Taylor instability with free surface"
set xlabel 'Time [My]' offset 0.5,0,0 font "Times-New-Roman:Bold, 12"
set ylabel 'Drip depth [km]' offset 0,0.5,0 font "Times-New-Roman:Bold, 12"
set output 'drip_depth.png'

set style line 1 lc rgb 'black' lt 1 lw 2 pt 1 ps 0.8
set style line 2 lc rgb 'brown' lt 1 lw 2 pt 2 pi -1 ps 0.8
set style line 3 lc rgb '#D2691E' lt 1 lw 2 pt 3 ps 0.8
set style line 4 lc rgb 'orange' lt 1 lw 2 pt 4 pi -1 ps 0.8
set style line 5 lc rgb 'black' lt 1 lw 2 pt 7 pi -1 ps 0.8

set xrange [0:6]
set yrange [480:80]

set key bottom right

plot 'output_nostab_2500/statistics' every 30 u ($2/1e6):($19/1e3) ls 1 w lp t 'no stab., dt=2500 yr', \
     'output_nostab_10000/statistics' every 35 u ($2/1e6):($19/1e3) ls 2 w lp t 'no stab., dt=10000 yr', \
     'output_stab_10000/statistics' every 40 u ($2/1e6):($19/1e3) ls 3 w lp t 'stab., dt=10000 yr', \
     'output_stab_10000_normal/statistics' every 45 u ($2/1e6):($19/1e3) ls 4 w lp t 'stab., dt=10000 yr, normal', \
     'Kaus_2010.csv' u ($1):(-$2) ls 5 w l t 'Kaus et al. 2010, dt=2500 yr'
