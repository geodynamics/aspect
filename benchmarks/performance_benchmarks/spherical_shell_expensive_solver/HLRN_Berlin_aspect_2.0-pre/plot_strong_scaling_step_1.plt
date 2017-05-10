set terminal pdf color dashed enhanced font 'Arial'

set output "plot.pdf"
 
set datafile missing '0'

set grid ytics

set xlabel "#Cores"
set ylabel "Wallclock time [s]"
set logscale xy

set xrange [1:15000]
set yrange [0.1:1000]

set key bottom left noautotitles

plot "runtimes_step_1.txt" using 1:2 with lines lt 1 lc rgb "blue" lw 3 title "Resolution 2, DoFs 2.7e5", \
     "runtimes_step_1.txt" using 1:3 with lines lt 1 lc rgb "orange" lw 3 title "Resolution 3, DoFs 2.1e6", \
     "runtimes_step_1.txt" using 1:4 with lines lt 1 lc rgb "yellow2" lw 3 title "Resolution 4, DoFs 1.6e7", \
     "runtimes_step_1.txt" using 1:5 with lines lt 1 lc rgb "green" lw 3 title "Resolution 5, DoFs 1.3e8", \
     "runtimes_step_1.txt" using 1:6 with lines lt 1 lc rgb "dark-red" lw 3 title "Resolution 6, DoFs 1.0e9", \
     1000/x lt 1 lc rgb "light-blue" lw 3 title "Optimal"
#     "runtimes_step_1_one_leaf.txt" using 1:2 with points lc rgb "blue" , \
#     "runtimes_step_1_one_leaf.txt" using 1:3 with points lc rgb "orange", \
#     "runtimes_step_1_one_leaf.txt" using 1:4 with points lc rgb "yellow2", \

