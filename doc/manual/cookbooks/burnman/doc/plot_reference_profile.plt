# Plot the reference state

set ylabel "Depth in km"
set xlabel "Temperature in K"
unset key 
set yrange [2800:0.0]
set datafile missing "nan"
set xtics (1600, 2000, 2400)

set terminal pdf color dashed enhanced font 'Arial,11' size 20cm,7.5cm
set output 'reference_profile.pdf'

set multiplot layout 1,6

set bmargin at screen 0.15
set lmargin 9
set rmargin 1

set title "Reference temperature"
plot "../../../../data/adiabatic-conditions/ascii-data/isentrope_properties.txt" using 3:($1/1000) with lines linecolor rgb "gray" lw 4 lt 1 title "temperature"
	
set lmargin 1

set xtics (3500, 4000, 4500, 5000)
set xlabel "Density in kg/m^3"
unset ylabel
set ytics format " " 

set title "Reference density"
plot "../../../../data/adiabatic-conditions/ascii-data/isentrope_properties.txt" using 4:($1/1000) with lines linecolor rgb "gray" lw 4 lt 1 title "density"

set xlabel "Gravity in m/s^2"
set style fill pattern 6 border
set xtics (9.75, 10.0, 10.25, 10.5)

set title "Gravity profile"
plot "../../../../data/adiabatic-conditions/ascii-data/isentrope_properties.txt" using 5:($1/1000) with lines linecolor rgb "gray" lw 4 lt 1 title "gravity"

set xlabel "Thermal expansivity in 1/K"
set xtics (2e-5, 4e-5, 6e-5)
set title "Thermal expansivity"
plot "../../../../data/adiabatic-conditions/ascii-data/isentrope_properties.txt" using 6:($1/1000) with lines linecolor rgb "gray" lw 4 lt 1 title "thermal expansivity"

set xlabel "Specific heat in J/kg/K"
set xtics (1220, 1240, 1260)
set title "Specific heat"
plot "../../../../data/adiabatic-conditions/ascii-data/isentrope_properties.txt" using 7:($1/1000) with lines linecolor rgb "gray" lw 4 lt 1 title "specific heat"

set xlabel "Compressibility in 1/Pa"
set xtics (3e-12, 6e-12, 9e-12, 1.2e-11)
set format x '%6.6g'
set title "Compressibility"
plot "../../../../data/adiabatic-conditions/ascii-data/isentrope_properties.txt" using 8:($1/1000) with lines linecolor rgb "gray" lw 4 lt 1 title "compressibility"

unset multiplot
