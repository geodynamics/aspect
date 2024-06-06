#### Script for plotting the results of the Quinquis 
#### subduction benchmark as described in the thesis
#### Quinquis, M. (2014). A numerical study of subduction
#### zone dynamics using linear viscous to thermo-mechanical
#### model setups including (de)hydration processes. 
#### Charles University, Prague.
#### Note that some of the results were provided by
#### Cedric Thieulot who ran the same benchmarks with
#### the finite element code Elefant.

set term png 
set output "Case2b_diagnostics.png"

#### Plot 2x3 graphs on one page 
set multiplot layout 3,3 rowsfirst

#### Some plotting settings that are the same for every graph
set xtics 2
set xtics font ", 10" 
set xlabel 'Time [My]' font ", 10"
set xrange [0:15]
set ytics font ", 10" 
set key font ",8"
set key samplen 2 
set key spacing 1
set grid lw 1
set size square

#### x- and ytics for each row resp. column
NOXTICS = "set xtics ('' 0, '' 2, '' 4,'' 6,'' 8,'' 10,'' 12,'' 14); unset xlabel"
XTICS = "set xtics 2; set xlabel 'Time [My]' font ', 10'"
#### Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.97; set bmargin at screen 0.73"
MMARGIN = "set tmargin at screen 0.67; set bmargin at screen 0.43"
BMARGIN = "set tmargin at screen 0.37; set bmargin at screen 0.13"
LMARGIN = "set lmargin at screen 0.13; set rmargin at screen 0.33"
CMARGIN = "set lmargin at screen 0.45; set rmargin at screen 0.65"
RMARGIN = "set lmargin at screen 0.78; set rmargin at screen 0.98"

#### The path to the ASPECT output
case = "output-Case2b/"
file = case . "statistics"


#### Linestyles
set style line 1 dt 1 lw 1.6 lc rgb "blue"

##### Slab tip depth ####
set key left bottom
@NOXTICS
@TMARGIN
@LMARGIN
set yrange [-670.:-80.]
set ytics 100 font ",10"
set ylabel 'Slab tip depth [km]' font ",10"
plot file u ($2/1e6):(-$27/1e3) w l ls 1 t 'ASPECT'
 

#### Trench location ####
@TMARGIN
@CMARGIN
unset key
set yrange [1340.:1480.]
set ytics 20 font ",10"
set ylabel 'Trench position [km]' font ",10"
plot file u ($2/1e6):($49/1e3) w l ls 1 notitle


#### Isotherm depth ####
@TMARGIN
@RMARGIN
set yrange [-670.:0.]
set ytics 100 font ",10"
set ylabel 'Isotherm depth [km]' font ",10"
plot file u ($2/1e6):(-$50/1e3) w l ls 1 notitle


#### Total Vrms ####
@MMARGIN
@LMARGIN
set yrange [2.0:5.0]
set ytics 0.5
set ylabel 'Vrms [cm/yr]' font ", 10"
plot file u ($2/1e6):($20*100.) w l ls 1 notitle
 

#### Total viscous dissipation ####
@MMARGIN
@CMARGIN
set yrange [6:20]
set ytics 2.0
set ylabel 'Total D [kW/m]' font ", 10"
plot file u ($2/1e6):($44/1e3) w l ls 1 notitle
 

#### Total Trms ####
@MMARGIN
@RMARGIN
set yrange [1286:1310]
set ytics 4.0
set ylabel 'Total Trms [C]' font ", 10"
plot file u ($2/1e6):($59) w l ls 1 notitle
 
#### Slab Vrms ####
@BMARGIN
@LMARGIN
@XTICS
set yrange [4.0:5.2]
set ytics 0.2
set ylabel 'Slab Vrms [cm/yr]' font ", 10"
plot file u ($2/1e6):($36*100.) w l ls 1 notitle

 
#### Slab viscous dissipation ####
@BMARGIN
@CMARGIN
set yrange [3.0:9.0]
set ytics 1.0
set ylabel 'Slab D [kW/m]' font ", 10"
plot file u ($2/1e6):(($38+$40+$42)/1e3) w l ls 1 notitle


#### Slab Trms ####
@BMARGIN
@RMARGIN
set yrange [820:855]
set ytics 5
set ylabel 'Slab Trms [C]' font ", 10"
plot file u ($2/1e6):($58) w l ls 1 notitle

unset multiplot
