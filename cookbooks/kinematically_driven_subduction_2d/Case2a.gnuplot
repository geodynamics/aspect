#### Script for plotting the results of the Quinquis 
#### subduction benchmark as described in the thesis
#### Quinquis, M. (2014). A numerical study of subduction
#### zone dynamics using linear viscous to thermo-mechanical
#### model setups including (de)hydration processes. 
#### Charles University, Prague.

set term png 
set output "Case2a_diagnostics_test2023.png"

#### Plot 3x3 graphs on one page 
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

#### The path to the directory
dir = "../../../HLRN/HLRN/QQ_subduction/"
case = "Case2/Case2a/"
file = dir . case . "statistics"

#### Data from other codes
SULEC = "SULEC_Case1_diagnostics.dat"
ELEFANT = "Elefant_Case2a_diagnostics.dat"
ELEFANTstd = "Elefant_Case2a_diagnostics_slabtipdepth.dat"
ELEFANTt = "Elefant_Case2a_diagnostics_trench.dat"

#### Linestyles
set style line 1 dt 1 lw 2 lc rgb "black"
set style line 2 dt 2 lw 1 lc rgb "black"

##### Slab tip depth ####
set key left bottom
@NOXTICS
@TMARGIN
@LMARGIN
set yrange [-670.:-80.]
set ytics 100 font ",10"
set ylabel 'Slab tip depth [km]' font ",10"
plot file u ($2/1e6):(-$22) ls 1 w l t 'ASPECT', \
     "Elefant_Case2a_tip.dat" u ($2/1.e6):(-670+$5/1000) ls 4 w l t 'Elefant'
 
##### Trench location ####
@TMARGIN
@CMARGIN
unset key
set yrange [1340.:1480.]
set ytics 20 font ",10"
set ylabel 'Trench position [km]' font ",10"
plot file u ($2/1e6):($23) w l ls 1 notitle, \
     "Elefant_Case2a_trench.dat" u ($2/1.e6):($4/1000) ls 4 w l notitle
 
##### Isotherm depth ####
@TMARGIN
@RMARGIN
set yrange [-670.:0.]
set ytics 100 font ",10"
set ylabel 'Isotherm depth [km]' font ",10"
plot file u ($2/1e6):(-$29) w l ls 1 notitle, \
     "Elefant_Case2a_isotherm.dat" u ($1/1.e6):($2) ls 4 w l notitle

#### Total Vrms ####
@MMARGIN
@LMARGIN
set yrange [2.0:5.0]
set ytics 0.5
set ylabel 'Vrms [cm/yr]' font ", 10"
plot file u ($2/1e6):($20*100.) w l ls 1 notitle, \
          "Elefant_Case2a_diagnostics_vrms.dat" u ($2/1e6):(100*$4) ls 4 w l notitle
 
#### Total viscous dissipation ####
#### TODO: add ASPECT dissipation
@MMARGIN
@CMARGIN
set yrange [6:20]
set ytics 2.0
set ylabel 'Total D [kW/m]' font ", 10"
plot "Elefant_Case2a_diagnostics_viscdiss.dat" u ($2/1e6):($4/1000) ls 4 w l notitle
#plot file u ($2/1e6):($28) w l ls 1 notitle, \
 
##### Total Trms ####
#### TODO: add ASPECT RMS temperature
@MMARGIN
@RMARGIN
set yrange [1287:1310]
set ytics 10.0
set ylabel 'Total Trms [C]' font ", 10"
plot "Elefant_Case2a_diagnostics_Trms.dat"  u ($2/1e6):4 ls 4 w l notitle
#plot file u ($2/1e6):($24) w l ls 1 notitle, \
 
##### Slab Vrms ####
@BMARGIN
@LMARGIN
@XTICS
set yrange [4.0:5.2]
set ytics 0.2
set ylabel 'Slab Vrms [cm/yr]' font ", 10"
plot file u ($2/1e6):($26) w l ls 1 notitle, \
          "Elefant_Case2a_diagnostics_vrms.dat" u ($2/1e6):(100*$3) ls 4 w l notitle
 
##### Slab viscous dissipation ####
@BMARGIN
@CMARGIN
set yrange [3.0:9.0]
set ytics 1.0
set ylabel 'Slab D [kW/m]' font ", 10"
plot file u ($2/1e6):($28) w l ls 1 notitle, \
     "Elefant_Case2a_diagnostics_viscdiss.dat"  u ($2/1e6):($3/1000) ls 4 w l notitle

##### Slab Trms ####
@BMARGIN
@RMARGIN
set yrange [816:850]
set ytics 5
set ylabel 'Slab Trms [C]' font ", 10"
plot file u ($2/1e6):($24) w l ls 1 notitle, \
     "Elefant_Case2a_diagnostics_Trms.dat"  u ($2/1e6):3 ls 4 w l notitle
 
unset multiplot
