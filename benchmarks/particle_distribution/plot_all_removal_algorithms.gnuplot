plot "output-cutoff-w1/statistics" using 2:23 title "deal.ii-cutoff-w1", \
     "output-cutoff-c1/statistics" using 2:23 title "deal.ii-cutoff-c1", \
     "output-gaussian/statistics" using 2:23 title "gaussian", \
     "output-triangular/statistics" using 2:23 title "triangular", \
     "output-uniform/statistics" using 2:23 title "uniform", \
     "output-random/statistics" using 2:23 title "random"
set style data linespoints
set xl "Time"
set ylabel "Standard Deviation of Cell Scores"

set terminal png
set output "cell-stdev-mean-oscillate.png"
replot



plot "output-cutoff-w1-constant-velocity/statistics" using 2:23 title "deal.ii-cutoff-w1 constant V", \
     "output-cutoff-c1-constant-velocity/statistics" using 2:23 title "deal.ii-cutoff-c1 constant V", \
     "output-gaussian-constant-velocity/statistics" using 2:23 title "gaussian constant V", \
     "output-triangular-constant-velocity/statistics" using 2:23 title "triangular constant V", \
     "output-uniform-constant-velocity/statistics" using 2:23 title "uniform constant V", \
     "output-random-constant-velocity/statistics" using 2:23 title "random constant V"
set style data linespoints
set xl "Time"
set ylabel "Standard Deviation of Cell Scores"

set terminal png
set output "cell-stdev-mean-constant.png"
replot



plot "output-cutoff-w1/statistics" using 2:21 title "deal.ii-cutoff-w1", \
     "output-cutoff-c1/statistics" using 2:21 title "deal.ii-cutoff-c1", \
     "output-gaussian/statistics" using 2:21 title "gaussian", \
     "output-triangular/statistics" using 2:21 title "triangular", \
     "output-uniform/statistics" using 2:21 title "uniform", \
     "output-random/statistics" using 2:21 title "random"
set style data linespoints
set xl "Time"
set ylabel "Average Cell Score"

set terminal png
set output "ave-score-oscillate.png"
replot



plot "output-cutoff-w1-constant-velocity/statistics" using 2:21 title "deal.ii-cutoff-w1 constant V", \
     "output-cutoff-c1-constant-velocity/statistics" using 2:21 title "deal.ii-cutoff-c1 constant V", \
     "output-gaussian-constant-velocity/statistics" using 2:21 title "gaussian constant V", \
     "output-triangular-constant-velocity/statistics" using 2:21 title "triangular constant V", \
     "output-uniform-constant-velocity/statistics" using 2:21 title "uniform constant V", \
     "output-random-constant-velocity/statistics" using 2:21 title "random constant V"
set style data linespoints
set xl "Time"
set ylabel "Average Cell Score"

set terminal png
set output "ave-score-constant.png"
replot



plot "output-cutoff-w1/statistics" using 2:22 title "deal.ii-cutoff-w1", \
     "output-cutoff-c1/statistics" using 2:22 title "deal.ii-cutoff-c1", \
     "output-gaussian/statistics" using 2:22 title "gaussian", \
     "output-triangular/statistics" using 2:22 title "triangular", \
     "output-uniform/statistics" using 2:22 title "uniform", \
     "output-random/statistics" using 2:22 title "random"
set style data linespoints
set xl "Time"
set ylabel "Maximum Cell Score"

set terminal png
set output "max-score-oscillate.png"
replot



plot "output-cutoff-w1-constant-velocity/statistics" using 2:22 title "deal.ii-cutoff-w1 constant V", \
     "output-cutoff-c1-constant-velocity/statistics" using 2:22 title "deal.ii-cutoff-c1 constant V", \
     "output-gaussian-constant-velocity/statistics" using 2:22 title "gaussian constant V", \
     "output-triangular-constant-velocity/statistics" using 2:22 title "triangular constant V", \
     "output-uniform-constant-velocity/statistics" using 2:22 title "uniform constant V", \
     "output-random-constant-velocity/statistics" using 2:22 title "random constant V"
set style data linespoints
set xl "Time"
set ylabel "Maximum Cell Score"

set terminal png
set output "max-score-constant.png"
replot
