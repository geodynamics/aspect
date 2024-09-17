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
set output "Case2a_Case1_diagnostics.png"

#### Plot 3x3 graphs on one page 
set multiplot layout 3,3 rowsfirst

#### Some plotting settings that are the same for every graph
set xtics 2
set xtics font ", 10" 
set xlabel 'Time [My]' font ", 10"
set xrange [0:15]
set ytics font ", 10" 
set key font ",7"
set key samplen 0.6 
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
case1 = "output-Case1/"
case2 = "output-Case2a/"
file1 = case1 . "statistics"
file2 = case2 . "statistics"

##### Data from other codes
ELEFANT = "Elefant/"

#### Linestyles
set style line 1 dt 1 lw 1 lc rgb "blue" pt 4 ps 0.8 pn 15
set style line 2 dt 1 lw 1.6 lc rgb "black"
set style line 3 dt 1 lw 1 lc rgb "orange"

##### Slab tip depth ####
@NOXTICS
@TMARGIN
@LMARGIN
set yrange [-670.:-80.]
set ytics 100 font ",10"
set ylabel 'Slab tip depth [km]' font ",10"
plot file1 u ($2/1e6):(-$27/1e3) w lp ls 1 notitle, \
     file2 u ($2/1e6):(-$27/1e3) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_tip.dat" u ($2/1.e6):(-670+$5/1000) w l ls 3 notitle
 
##### Trench location ####
@TMARGIN
@CMARGIN
unset key
set yrange [1390.:1480.]
set ytics 10 font ",10"
set ylabel 'Trench position [km]' font ",10"
plot file1 u ($2/1e6):($45/1e3) w lp ls 1 notitle, \
     file2 u ($2/1e6):($49/1e3) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_trench.dat" u ($2/1.e6):($4/1000) w l ls 3 notitle
 
###### Isotherm depth ####
@TMARGIN
@RMARGIN
set yrange [-670.:0.]
set ytics 100 font ",10"
set ylabel 'Isotherm depth [km]' font ",10"
plot file2 u ($2/1e6):(-$50/1000.) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_isotherm.dat" u ($1/1.e6):($2) w l ls 3 notitle

#### Total Vrms ####
@MMARGIN
@LMARGIN
set yrange [2.0:3.4]
set ytics 0.2
set ylabel 'Vrms [cm/yr]' font ", 10"
plot file1 u ($2/1e6):($20*100.) w lp ls 1 notitle, \
     file2 u ($2/1e6):($20*100.) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_diagnostics_vrms.dat" u ($2/1e6):(100*$4) w l ls 3 notitle
 
#### Total viscous dissipation ####
set key left top spacing 1.2 font "Arial, 7"
@MMARGIN
@CMARGIN
set yrange [6:26]
set ytics 2.0
set ylabel 'Total D [kW/m]' font ", 10"
plot file1 u ($2/1e6):($44/1e3) w lp ls 1 t 'ASPECT-case1', \
     file2 u ($2/1e6):($44/1e3) w l ls 2 t 'ASPECT-case2a', \
     ELEFANT . "Elefant_Case2a_diagnostics_viscdiss.dat" u ($2/1e6):($4/1000) w l ls 3 t 'Elefant-case2a'
 
###### Total Trms ####
@MMARGIN
@RMARGIN
set yrange [1286:1310]
set ytics 4.0
set ylabel 'Total Trms [C]' font ", 10"
plot file2 u ($2/1e6):($59) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_diagnostics_Trms.dat" u ($2/1e6):4 w l ls 3 notitle
 
##### Slab Vrms ####
@BMARGIN
@LMARGIN
@XTICS
set yrange [3.6:6.0]
set ytics 0.4
set ylabel 'Slab Vrms [cm/yr]' font ", 10"
plot file1 u ($2/1e6):($36*100.) w lp ls 1 notitle, \
     file2 u ($2/1e6):($36*100.) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_diagnostics_vrms.dat" u ($2/1e6):(100*$3) w l ls 3 notitle
 
##### Slab viscous dissipation ####
@BMARGIN
@CMARGIN
set yrange [2.0:18.0]
set ytics 2.0
set ylabel 'Slab D [kW/m]' font ", 10"
plot file1 u ($2/1e6):(($38+$40+$42)/1e3) w lp ls 1 notitle, \
     file2 u ($2/1e6):(($38+$40+$42)/1e3) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_diagnostics_viscdiss.dat" u ($2/1e6):($3/1000) w l ls 3 notitle

###### Slab Trms ####
@BMARGIN
@RMARGIN
set yrange [820:870]
set ytics 5
set ylabel 'Slab Trms [C]' font ", 10"
plot file2 u ($2/1e6):($58) w l ls 2 notitle, \
     ELEFANT . "Elefant_Case2a_diagnostics_Trms.dat"  u ($2/1e6):3 w l ls 3 notitle
 
unset multiplot
