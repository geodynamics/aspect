plot "output_addition/output-cutoff-w1/statistics" using 2:23 title "deal.ii-cutoff-w1", \
     "output_addition/output-cutoff-c1/statistics" using 2:23 title "deal.ii-cutoff-c1", \
     "output_addition/output-gaussian/statistics" using 2:23 title "gaussian", \
     "output_addition/output-triangular/statistics" using 2:23 title "triangular", \
     "output_addition/output-uniform/statistics" using 2:23 title "uniform", \
     "output_addition/output-histogram/statistics" using 2:23 title "histogram", \
     "output_addition/output-random/statistics" using 2:23 title "random"
set style data linespoints
set xl "Time"
set ylabel "Standard Deviation of Cell Scores"

set terminal png
set output "./doc/cell-stdev-mean-oscillate-addition.png"
replot



plot "output_addition/output-cutoff-w1-constant-velocity/statistics" using 2:23 title "deal.ii-cutoff-w1 constant V", \
     "output_addition/output-cutoff-c1-constant-velocity/statistics" using 2:23 title "deal.ii-cutoff-c1 constant V", \
     "output_addition/output-gaussian-constant-velocity/statistics" using 2:23 title "gaussian constant V", \
     "output_addition/output-triangular-constant-velocity/statistics" using 2:23 title "triangular constant V", \
     "output_addition/output-uniform-constant-velocity/statistics" using 2:23 title "uniform constant V", \
     "output_addition/output-histogram-constant-velocity/statistics" using 2:23 title "histogram constant V", \
     "output_addition/output-random-constant-velocity/statistics" using 2:23 title "random constant V"
set style data linespoints
set xl "Time"
set ylabel "Standard Deviation of Cell Scores"

set terminal png
set output "./doc/cell-stdev-mean-constant-addition.png"
replot



plot "output_addition/output-cutoff-w1/statistics" using 2:21 title "deal.ii-cutoff-w1", \
     "output_addition/output-cutoff-c1/statistics" using 2:21 title "deal.ii-cutoff-c1", \
     "output_addition/output-gaussian/statistics" using 2:21 title "gaussian", \
     "output_addition/output-triangular/statistics" using 2:21 title "triangular", \
     "output_addition/output-uniform/statistics" using 2:21 title "uniform", \
     "output_addition/output-histogram/statistics" using 2:21 title "histogram", \
     "output_addition/output-random/statistics" using 2:21 title "random"
set style data linespoints
set xl "Time"
set ylabel "Average Cell Score"

set terminal png
set output "./doc/ave-score-oscillate-addition.png"
replot



plot "output_addition/output-cutoff-w1-constant-velocity/statistics" using 2:21 title "deal.ii-cutoff-w1 constant V", \
     "output_addition/output-cutoff-c1-constant-velocity/statistics" using 2:21 title "deal.ii-cutoff-c1 constant V", \
     "output_addition/output-gaussian-constant-velocity/statistics" using 2:21 title "gaussian constant V", \
     "output_addition/output-triangular-constant-velocity/statistics" using 2:21 title "triangular constant V", \
     "output_addition/output-uniform-constant-velocity/statistics" using 2:21 title "uniform constant V", \
     "output_addition/output-histogram-constant-velocity/statistics" using 2:21 title "histogram", \
     "output_addition/output-random-constant-velocity/statistics" using 2:21 title "random constant V"
set style data linespoints
set xl "Time"
set ylabel "Average Cell Score"

set terminal png
set output "./doc/ave-score-constant-addition.png"
replot



plot "output_addition/output-cutoff-w1/statistics" using 2:22 title "deal.ii-cutoff-w1", \
     "output_addition/output-cutoff-c1/statistics" using 2:22 title "deal.ii-cutoff-c1", \
     "output_addition/output-gaussian/statistics" using 2:22 title "gaussian", \
     "output_addition/output-triangular/statistics" using 2:22 title "triangular", \
     "output_addition/output-uniform/statistics" using 2:22 title "uniform", \
     "output_addition/output-histogram/statistics" using 2:22 title "histogram", \
     "output_addition/output-random/statistics" using 2:22 title "random"
set style data linespoints
set xl "Time"
set ylabel "Maximum Cell Score"

set terminal png
set output "./doc/max-score-oscillate-addition.png"
replot



plot "output_addition/output-cutoff-w1-constant-velocity/statistics" using 2:22 title "deal.ii-cutoff-w1 constant V", \
     "output_addition/output-cutoff-c1-constant-velocity/statistics" using 2:22 title "deal.ii-cutoff-c1 constant V", \
     "output_addition/output-gaussian-constant-velocity/statistics" using 2:22 title "gaussian constant V", \
     "output_addition/output-triangular-constant-velocity/statistics" using 2:22 title "triangular constant V", \
     "output_addition/output-uniform-constant-velocity/statistics" using 2:22 title "uniform constant V", \
     "output_addition/output-histogram-constant-velocity/statistics" using 2:22 title "histogram", \
     "output_addition/output-random-constant-velocity/statistics" using 2:22 title "random constant V"
set style data linespoints
set xl "Time"
set ylabel "Maximum Cell Score"

set terminal png
set output "./doc/max-score-constant-addition.png"
replot
