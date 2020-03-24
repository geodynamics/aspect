# Plot the root mean square velocity

set ylabel "Root mean square velocity in m/yr"
set xlabel "Times in years"
set key bottom 
set yrange [0.015:0.04]
set xrange [0.0:2.63e+08]
set datafile missing "nan"

set terminal pdf color dashed enhanced font 'Arial,11' size 9cm,7.5cm
set output 'vrms.pdf'

set title "Comparison between the different approximations"
plot "statistics_ala" using 2:12 with lines lw 3 linecolor rgb "navy" title "Anelastic liquid approximation", \
     "statistics_tala" using 2:12 with lines lt 1 lw 3 linecolor rgb "web-green" title "Truncated anelastic liquid approximation", \
     "statistics_ica" using 2:12 with lines lt 1 lw 3 linecolor rgb "light-red" title "Isothermal compression approximation"



