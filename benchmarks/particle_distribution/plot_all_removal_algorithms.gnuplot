plot "output-cutoff-w1/statistics" using 2:26 title "deal.ii-cutoff-w1", \
     "output-gaussian/statistics" using 2:26 title "gaussian", \
     "output-triangular/statistics" using 2:26 title "triangular", \
     "output-uniform/statistics" using 2:26 title "uniform", \
     "output-random/statistics" using 2:26 title "random"
set style data linespoints
set xl "Time"
set ylabel "stdev mean"

set terminal png
set output "stdev-mean-oscillate.png"
replot



plot "output-cutoff-w1-constant-velocity/statistics" using 2:26 title "deal.ii-cutoff-w1 constant V", \
     "output-gaussian-constant-velocity/statistics" using 2:26 title "gaussian constant V", \
     "output-triangular-constant-velocity/statistics" using 2:26 title "triangular constant V", \
     "output-uniform-constant-velocity/statistics" using 2:26 title "uniform constant V", \
     "output-random-constant-velocity/statistics" using 2:26 title "random constant V"
set style data linespoints
set xl "Time"
set ylabel "stdev mean"

set terminal png
set output "stdev-mean-constant.png"
replot



plot "output-cutoff-w1/statistics" using 2:27 title "deal.ii-cutoff-w1", \
     "output-gaussian/statistics" using 2:27 title "gaussian", \
     "output-triangular/statistics" using 2:27 title "triangular", \
     "output-uniform/statistics" using 2:27 title "uniform", \
     "output-random/statistics" using 2:27 title "random" 
     
set style data linespoints
set xl "Time"
set ylabel "stdev max"

set terminal png
set output "stdev-max-oscillate.png"
replot


plot "output-cutoff-w1-constant-velocity/statistics" using 2:26 title "deal.ii-cutoff-w1 constant V", \
     "output-gaussian-constant-velocity/statistics" using 2:26 title "gaussian constant V", \
     "output-triangular-constant-velocity/statistics" using 2:26 title "triangular constant V", \
     "output-uniform-constant-velocity/statistics" using 2:26 title "uniform constant V", \
     "output-random-constant-velocity/statistics" using 2:26 title "random constant V"

set style data linespoints
set xl "Time"
set ylabel "stdev max"

set terminal png
set output "stdev-max-constant.png"
replot


plot "output-cutoff-w1/statistics" using 2:23 title "deal.ii-cutoff-w1", \
     "output-gaussian/statistics" using 2:23 title "gaussian", \
     "output-triangular/statistics" using 2:23 title "triangular", \
     "output-uniform/statistics" using 2:23 title "uniform", \
     "output-random/statistics" using 2:23 title "random"

     
set style data linespoints
set xl "Time"
set ylabel "score mean"

set terminal png
set output "score-mean-oscillate.png"
replot



plot "output-cutoff-w1-constant-velocity/statistics" using 2:26 title "deal.ii-cutoff-w1 constant V", \
     "output-gaussian-constant-velocity/statistics" using 2:26 title "gaussian constant V", \
     "output-triangular-constant-velocity/statistics" using 2:26 title "triangular constant V", \
     "output-uniform-constant-velocity/statistics" using 2:26 title "uniform constant V", \
     "output-random-constant-velocity/statistics" using 2:26 title "random constant V"

set style data linespoints
set xl "Time"
set ylabel "score mean"

set terminal png
set output "score-mean-constant.png"
replot



plot "output-cutoff-w1/statistics" using 2:24 title "deal.ii-cutoff-w1", \
     "output-gaussian/statistics" using 2:24 title "gaussian", \
     "output-triangular/statistics" using 2:24 title "triangular", \
     "output-uniform/statistics" using 2:24 title "uniform", \
     "output-random/statistics" using 2:24 title "random" 
     
set style data linespoints
set xl "Time"
set ylabel "score max"

set terminal png
set output "score-max-oscillate.png"
replot


plot "output-cutoff-w1-constant-velocity/statistics" using 2:26 title "deal.ii-cutoff-w1 constant V", \
     "output-gaussian-constant-velocity/statistics" using 2:26 title "gaussian constant V", \
     "output-triangular-constant-velocity/statistics" using 2:26 title "triangular constant V", \
     "output-uniform-constant-velocity/statistics" using 2:26 title "uniform constant V", \
     "output-random-constant-velocity/statistics" using 2:26 title "random constant V"

set style data linespoints
set xl "Time"
set ylabel "score max"

set terminal png
set output "score-max-constant.png"
replot
