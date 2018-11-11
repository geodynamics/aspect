!/bin/bash
declare -a version="c3"
#declare -a boundary_type=("velocities") # "tractions")
#for i_boundary_type in "${boundary_type[@]}"
#do
input_file="input_""$i_boundary_type"".prm"
declare -a grid=(4) #4 5 6 7 #2 3 4 5)
agrid=0
declare -a NSP=("SPD") # "symmetric" "PD" "none") #Use Newton stabilization preconditioner
declare -a NSA=("SPD") # "symmetric" "none") # "Symmetric" "PD" "none") #Use Newton stabilization A block
declare -a UFS=("false" "true")
declare -a n=(50) # 5) #3 5 #1 2 5 10 25 50 100)
I=150
declare -a P=(0 5 10 15 20 25 30 50 150) # 5 10 100) #0 5 10 100 
declare -a ST=(1e-20)
declare -a LS=(0 5 10) # 5 50) # 5 10 50) # 5 10) # 0,5,10
declare -a LT=("1e-5") # "1e-8") # 6 7 8) # 5,6,7,8
declare -a NLT=("1e-14")
declare -a RSM=("false" "true")
declare -a UDS=("true") #"true" #"false"
declare -a SF=("9e-1" "9.9e-1") # "9.99e-1" "1")
declare -a AEW=("false")
AV=-1
ATB="1e-2"
COMP="0" #"4e-12"
materialmodelnameShort="DP" #"VP2"
materialmodelname="DP"
phi=30
declare -a vel=(25 50 125)
declare -a BV=("1e23" "1e24" "5e24")
SOLVER_SHORT="NS" #"itIMPES"
SOLVER="iterated IMPES"
if [ $SOLVER_SHORT == "NS" ]; then
SOLVER="Newton Stokes"
elif [ $SOLVER_SHORT == "itIMPES" ]; then
SOLVER="iterated IMPES"
fi  
declare -a OS=("9e-1" "1e-1" "1e-2" "1e-8") #"1e-7" "1e-6" "1e-5" "1e-4" "1e-3" "1e-2" "1e-1" "5e-1" "9e-1") #"9e-1" "5e-1" "1e-1" "1e-2" "1e-4")

for i_UFS in "${UFS[@]}"
do
for i_NSP in "${NSP[@]}"
do
for i_NSA in "${NSA[@]}"
do
for i_ST in "${ST[@]}"
do
for i_RSM in "${RSM[@]}"
do
for i_grid in "${grid[@]}"
do
for i_LT in "${LT[@]}"
do
for i_LS in "${LS[@]}"
do
for i_UDS in "${UDS[@]}"
do
for i_SF in "${SF[@]}"
do
for i_AEW in "${AEW[@]}"
do
for i_NLT in "${NLT[@]}"
do
for i_n in "${n[@]}"
do

basename="$version""P-minLT-res-it-vel-BV_g""$i_grid""_ag""$agrid""_ST""$i_ST""_UFS""$i_UFS""_NSP-""$i_NSP""_NSA-""$i_NSA""_ATB""$ATB""_n""$i_n""_I""$I""_LS""$i_LS""_LT""$i_LT""_NLT""$i_NLT""_RSM""$i_RSM""_""$SOLVER_SHORT""_AV""$AV""_phi""$phi""_comp""$COMP""_UDS""$i_UDS""_SF""$i_SF""_AEW""$i_AEW"
#final_plot_file_name=${final_file_name::-4}".gnuplot"
final_plot_file_name="image_plotfiles/""$basename"".gnuplot"
final_png_file_name="image_collections/$basename.png"
rm $final_plot_file_name

echo "set terminal pngcairo size 3000,2400 enhanced font 'Verdana,20'" >> $final_plot_file_name
echo "set output '$final_png_file_name'" >> $final_plot_file_name
echo "set multiplot title \"{/*2.0 Final nonlinear iteration and final nonlinear residual as function of grid resolution and Max nonlinar tolerance}\\n\\n{/*1.5 Grid: $i_grid, LS: $i_LS, ST: $i_ST, NLT: $i_NLT, RSM: $i_RSM, phi: $phi, comp: $COMP, solver: $SOLVER_SHORT, UFS: $i_UFS, NSP: $i_NSP, NSA: $i_NSA, UDS: $i_UDS, SF: $i_SF, AEW: $i_AEW}\"" >> $final_plot_file_name
echo "set rmargin 2" >> $final_plot_file_name
echo "set bmargin 2.2" >> $final_plot_file_name
for i_BV in "${BV[@]}"
do
for i_vel in "${vel[@]}"
do

final_file_name="image_data/""$basename""_vel""$i_vel""_BV""$i_BV"".dat" #P-NLT-res-it_g""$i_grid""_n""$i_n""_I""$I""_LS""$i_LS""_LT""$i_LT""_minLT""$i_OS""_""$SOLVER_SHORT"".dat"
rm $final_file_name


for i_P in "${P[@]}"
do
for i_OS in "${OS[@]}"
do
#dirname="b3_""$i_boundary_type""_g""$i_grid""_n""$i_n""_I""$I""_P""$i_P""_LS""$i_LS""_LT""$i_LT""_NLT""$i_NLT""_UA""$i_UA""_minLT""$i_OS""_""$SOLVER_SHORT""_phi""$phi""_vel""$i_vel""_BV""i_BV"
dirname="results/""$version""$materialmodelnameShort""_""$SOLVER_SHORT""_ST""$i_ST""_UFS""$i_UFS""_NSP-""$i_NSP""_NSA-""$i_NSA""_C""$COMP""_g""$i_grid""_ag""$agrid""_AEW""$i_AEW""_UDS""$i_UDS""_SF""$i_SF""_NLT""$i_NLT""_ABT""$ATB""_LT""$i_LT""_mLT""$i_OS""_I""$I""_P""$i_P""_EW1""_theta1""_LS""$i_LS""_RSM""$i_RSM""_AV""$AV""_phi""$phi""_vel""$i_vel""_BV""$i_BV""_n""$i_n"

infilename="$dirname""/input.prm"
outfilename="$dirname""/output.log"
outplotfilename="$dirname""/plot.dat"
#infilename="g""$i_grid""_n""$i_n""_""$I""I_""$i_P""P_LS""$LS""_LT""$LT""/input.prm"
#outfilename="g""$i_grid""_n""$i_n""_""$I""I_""$i_P""P_LS""$LS""_LT""$LT""/output.log"
#outplotfilename="g""$i_grid""_n""$i_n""_""$I""I_""$i_P""P_LS""$LS""_LT""$LT""/plot.dat"
#dirname="g""${grid[i_grid]}""_n""${n[i_n]}""_""$I""I_""${P[i_P]}""P_LS""$LS""_LT""$LT"

#mkdir -p $dirname

echo "$dirname"

awk -v dirname="$dirname" -v var1="$i_P" -v var2="$i_OS" 'END{print dirname ", " var1 ", " var2 ", " substr($11,1,length($11)-1) ", " substr($10,1,length($10)-1)}' $outplotfilename >> $final_file_name
#sed  \
#-e "s/set Nonlinear solver scheme.*/set Nonlinear solver scheme = $SOLVER/g" \
#$input_file > "$infilename"

#mpirun -np 1 ../../.././aspect $infilename > $outfilename
#../../.././aspect $infilename > $outfilename

#grep "Total relative residual after nonlinear iteration " $outfilename > $outplotfilename
done
done
if [ "$i_vel" == "25" ] && [ "$i_BV" == "1e23" ]; then
   echo "set origin 0,0.64" >> $final_plot_file_name
   echo "set colorbox vert user origin .93,.0275 size .02,.92" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 25, BV = 1e23\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "50" ] && [ "$i_BV" == "1e23" ]; then
   echo "set origin 0.31,0.64" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 50, BV = 1e23\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "125" ] && [ "$i_BV" == "1e23" ]; then
   echo "set origin 0.62,0.64" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 125, BV = 1e23\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "25" ] && [ "$i_BV" == "1e24" ]; then
   echo "set origin 0,0.32" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 25, BV = 1e24\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "50" ] && [ "$i_BV" == "1e24" ]; then
   echo "set origin 0.31,0.32" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 50, BV = 1e24\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "125" ] && [ "$i_BV" == "1e24" ]; then
   echo "set origin 0.62,0.32" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 125, BV = 1e24\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "25" ] && [ "$i_BV" == "5e24" ]; then
   echo "set origin 0,0" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 25, BV = 5e24\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "50" ] && [ "$i_BV" == "5e24" ]; then
   echo "set origin 0.31,0" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 50, BV = 5e24\" offset 0,0.75" >> $final_plot_file_name
elif [ "$i_vel" == "125" ] && [ "$i_BV" == "5e24" ]; then
   echo "set origin 0.62,0" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
   echo "set xlabel \"Maximum Picard iterations, vel = 125, BV = 5e24\" offset 0,0.75" >> $final_plot_file_name
fi
echo "set size 0.29,0.32" >> $final_plot_file_name
echo "" >> $final_plot_file_name
echo "set ytics autofreq" >> $final_plot_file_name
echo "set ylabel \"Maximum nonlinear tolerance\" offset 1,0" >> $final_plot_file_name
echo "set cblabel \"Final nonlinear residual (square)\"" >> $final_plot_file_name
echo "set palette rgb 33,13,10" >> $final_plot_file_name
echo "set logscale y" >> $final_plot_file_name
echo "set logscale cb" >> $final_plot_file_name
echo "set format cb '%.0e'" >> $final_plot_file_name
echo "set cbtics nomirror" >> $final_plot_file_name
echo "" >> $final_plot_file_name
echo "set format y '%.0e'" >> $final_plot_file_name
echo "set xtics 10" >> $final_plot_file_name
echo "set mxtics 2" >> $final_plot_file_name
echo "set xrange [-3:53]" >> $final_plot_file_name
echo "set yrange [1e-9:10]" >> $final_plot_file_name
echo "set cbrange [1e-20:1e-1]" >> $final_plot_file_name
echo "unset key" >> $final_plot_file_name
echo "plot \"""$final_file_name""\" u 2:3:(3):4 w p ls 5 ps 10 lc palette " >> $final_plot_file_name
#echo "set colorbox horiz user origin .1,.07 size .7,.04" >> $final_plot_file_name
#echo "set cblabel \"Final nonlinear iteration (circle)\" offset 0,0.75" >> $final_plot_file_name
echo "unset colorbox" >> $final_plot_file_name
echo "unset logscale cb" >> $final_plot_file_name
echo "set format cb '%.0f'" >> $final_plot_file_name
echo "set cbtics nomirror" >> $final_plot_file_name
echo "set style fill solid" >> $final_plot_file_name
echo "set palette rgb 33,13,10" >> $final_plot_file_name
echo "set cbrange [0:100]" >> $final_plot_file_name
echo "set cbtics 20" >> $final_plot_file_name
echo "set mcbtics 4" >> $final_plot_file_name
#echo "plot \"""$final_file_name""\" u 2:3:(3):5 w p ls 7 ps 6 lc palette " >> $final_plot_file_name
#echo "plot \"""$final_file_name""\" u 2:3:(3):5 w p ls 6 ps 6 lc rgb \"black\" " >> $final_plot_file_name
echo "plot \"""$final_file_name""\" u 2:3:(sprintf('%.0f',\$5)) w labels font 'Times Bold,24' " >> $final_plot_file_name
#echo "plot \"""$final_file_name""\" u 2:3:(2) w circle lc rgb \"black\" " >> $final_plot_file_name
#echo "plot \"""$final_file_name""\" u 2:3:(1.5):5 w circle lc palette " >> $final_plot_file_name
#echo "pause -1" >> $final_plot_file_name

if [ "$i_vel" == "25" ] && [ "$i_BV" == "1e23" ]; then
   echo "set origin 0.265,0.64" >> $final_plot_file_name
elif [ "$i_vel" == "50" ] && [ "$i_BV" == "1e23" ]; then
   echo "set origin 0.575,0.64" >> $final_plot_file_name
elif [ "$i_vel" == "125" ] && [ "$i_BV" == "1e23" ]; then
   echo "set origin 0.885,0.64" >> $final_plot_file_name
   echo "unset colorbox" >> $final_plot_file_name
elif [ "$i_vel" == "25" ] && [ "$i_BV" == "1e24" ]; then
   echo "set origin 0.265,0.32" >> $final_plot_file_name
elif [ "$i_vel" == "50" ] && [ "$i_BV" == "1e24" ]; then
   echo "set origin 0.575,0.32" >> $final_plot_file_name
elif [ "$i_vel" == "125" ] && [ "$i_BV" == "1e24" ]; then
   echo "set origin 0.885,0.32" >> $final_plot_file_name
elif [ "$i_vel" == "25" ] && [ "$i_BV" == "5e24" ]; then
   echo "set origin 0.265,0" >> $final_plot_file_name
elif [ "$i_vel" == "50" ] && [ "$i_BV" == "5e24" ]; then
   echo "set origin 0.575,0" >> $final_plot_file_name
elif [ "$i_vel" == "125" ] && [ "$i_BV" == "5e24" ]; then
   echo "set origin 0.885,0" >> $final_plot_file_name
fi

echo "set size 0.05,0.32" >> $final_plot_file_name
echo "" >> $final_plot_file_name
echo "unset ylabel" >> $final_plot_file_name
echo "unset xlabel" >> $final_plot_file_name
echo "set cblabel \"Final nonlinear residual (square)\"" >> $final_plot_file_name
echo "set palette rgb 33,13,10" >> $final_plot_file_name
echo "set logscale y" >> $final_plot_file_name
echo "set logscale cb" >> $final_plot_file_name
echo "set format cb '%.0e'" >> $final_plot_file_name
echo "set cbtics nomirror" >> $final_plot_file_name
echo "" >> $final_plot_file_name
echo "set format y '%.0e'" >> $final_plot_file_name
echo "set xtics 10" >> $final_plot_file_name
echo "set mxtics 2" >> $final_plot_file_name
echo "unset ytics" >> $final_plot_file_name
echo "set xrange [147:153]" >> $final_plot_file_name
echo "set yrange [1e-9:10]" >> $final_plot_file_name
echo "set cbrange [1e-20:1e-1]" >> $final_plot_file_name
echo "unset key" >> $final_plot_file_name
echo "plot \"""$final_file_name""\" u 2:3:(3):4 w p ls 5 ps 10 lc palette " >> $final_plot_file_name
echo "unset colorbox" >> $final_plot_file_name
echo "unset logscale cb" >> $final_plot_file_name
echo "set format cb '%.0f'" >> $final_plot_file_name
echo "set cbtics nomirror" >> $final_plot_file_name
echo "set style fill solid" >> $final_plot_file_name
echo "set palette rgb 33,13,10" >> $final_plot_file_name
echo "set cbrange [0:100]" >> $final_plot_file_name
echo "set cbtics 20" >> $final_plot_file_name
echo "set mcbtics 4" >> $final_plot_file_name
echo "plot \"""$final_file_name""\" u 2:3:(sprintf('%.0f',\$5)) w labels font 'Times Bold,24' " >> $final_plot_file_name

gnuplot $final_plot_file_name
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
