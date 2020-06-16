#!/bin/bash
version="test"
build_dir="../../../build/"
declare -a grid=(5) 
declare -a UFS=("false" "true") # use failsafe
declare -a NSP=("SPD") # "Symmetric" "PD" "none") #Use Newton stabilization preconditioner
declare -a NSA=("SPD" "symmetric" "none") # "Symmetric" "PD" "none") #Use Newton stabilization A block
declare -a agrid=(0)
declare -a n=(50) #1 2 5 10 25 50 100)
I=150
declare -a P=(0 5 10 50 15 150 20 25 30) #0 5 10 15 20 150) # 1 2 3 4 5 10 15 25 50 100)
declare -a LS=(50)
declare -a ST=(1e-20)
declare -a LT=("1e-5") # "1e-8") # "1e-8") #"1e-5") # "1e-8")
declare -a NLT=("1e-14") # "1e-9" "1e-10" "1e-14" "1e-20")
declare -a ABT=("1e-2")
declare -a RSM=("true") #"true" #"false"
declare -a SF=("9e-1")
declare -a UDS=("true")
declare -a AEW=("false") # always use Eisenstat Walker, even for picard
AV=-1
COMP="0" #"4e-12"
declare -a OS=("9e-1" "1e-1" "1e-2" "1e-8") # "1e-3" "1e-4" "1e-5" "1e-6" "1e-6" "1e-7" "1e-8")
materialmodelnameShort="DP" #"VP2"
materialmodelname="DP"
phi=30 # 0
declare -a vel=(25 50 125) # 50 125) #25 50 125)
declare -a BV=("1e23" "1e24" "5e24") #"1e23" "1e24" "5e24")
processes=4

CohesionLine="s/set List of cohesions.*/      set List of cohesions                              = 1e8,0/g"
PhiLine="s/set List of angles of internal friction.*/      set List of angles of internal friction            = $phi,0/g"
if [ $materialmodelnameShort == "DP" ]; then
 materialmodelname="drucker prager compositions"
elif [ $materialmodelnameShort == "vM" ]; then
 materialmodelname="drucker prager compositions"
elif [ $materialmodelnameShort == "SNL" ]; then
 materialmodelname="simple nonlinear compositions"
 CohesionLine="s/ set List of cohesion of fields.*/ /g"
 PhiLine="s/set List of angles of internal friction of fields.*/ /g"
elif [ $materialmodelnameShort == "VP2" ]; then
 materialmodelname="viscoplastic2"
fi

SOLVER_SHORT="NS" #"itIMPES" #"NS" #NS"
SOLVER="NS"
if [ $SOLVER_SHORT == "NS" ]; then
 SOLVER="Newton Stokes"
elif [ $SOLVER_SHORT == "itIMPES" ]; then
 SOLVER="iterated IMPES"
fi


for i_grid in "${grid[@]}"
do
 for i_agrid in "${agrid[@]}"
 do
  for i_NSP in "${NSP[@]}"
  do
   for i_NSA in "${NSA[@]}"
   do
    for i_NLT in "${NLT[@]}"
    do
     for i_ABT in "${ABT[@]}"
     do
      for i_ST in "${ST[@]}"
      do
       for i_UDS in "${UDS[@]}"
       do
        for i_SF in "${SF[@]}"
        do
         for i_LT in "${LT[@]}"
         do
          for i_n in "${n[@]}"
          do
           for i_UFS in "${UFS[@]}"
           do
            for i_P in "${P[@]}"
            do
             for i_LS in "${LS[@]}"
             do
              for i_AEW in "${AEW[@]}"
              do
               for i_OS in "${OS[@]}"
               do
                for i_RSM in "${RSM[@]}"
                do
                 for i_BV in "${BV[@]}"
                 do
                  for i_vel in "${vel[@]}"
                  do
                  U="7.92219116e-11"
                  if [ $i_vel == 25 ]; then
                   U="7.92219116e-11"
                  elif [ $i_vel == 50 ]; then
                   U="1.58443823e-10"
                  elif [ $i_vel == 125 ]; then
                   U="3.96109558e-10"
                  fi
                  
                  dirname_clean="$version""$materialmodelnameShort""_""$SOLVER_SHORT""_ST""$i_ST""_UFS""$i_UFS""_NSP-""$i_NSP""_NSA-""$i_NSA""_C""$COMP""_g""$i_grid""_ag""$i_agrid""_AEW""$i_AEW""_UDS""$i_UDS""_SF""$i_SF""_NLT""$i_NLT""_ABT""$i_ABT""_LT""$i_LT""_mLT""$i_OS""_I""$I""_P""$i_P""_EW1""_theta1""_LS""$i_LS""_RSM""$i_RSM""_AV""$AV""_phi""$phi""_vel""$i_vel""_BV""$i_BV""_n""$i_n"
                  dirname="results/""$dirname_clean"
                  infilename="$dirname""/input.prm"
                  outfilename="$dirname""/output.log"
                  errorfilename="$dirname""/error.log"
                  outplotfilename="$dirname""/plot.dat"
                  
                  mkdir -p $dirname
                  
                  echo "$dirname"
                  sed  \
                   -e "$PhiLine" \
                   -e "$CohesionLine" \
                   -e "s/    set Function expression = if(x<60e3,.*/    set Function expression = if(x<60e3,$U,-$U);0/g" \
                   -e "s/    set Nonlinear Newton solver switch tolerance.*/     set Nonlinear Newton solver switch tolerance = $i_ST/g" \
                   -e "s/set Reference viscosity =.*/    set Reference viscosity = $i_BV/g" \
                   -e "s/set Output directory .*/set Output directory = results\/$dirname_clean/g" \
                   -e "s/    set Model name = .*/    set Model name = $materialmodelname/g" \
                   -e "s/set Nonlinear solver scheme.*/set Nonlinear solver scheme = $SOLVER/g" \
                   -e "s/  set Initial global refinement          = .*/  set Initial global refinement          = $i_grid/g" \
                   -e "s/  set Initial adaptive refinement        = .*/  set Initial adaptive refinement        = $i_agrid/g" \
                   -e "s/      set List of stress exponents of fields                       = .*/      set List of stress exponents of fields                       = $i_n, 1/g" \
                   -e "s/    set Viscosity averaging p = .*/    set Viscosity averaging p = $AV/g" \
                   -e "s/set Max nonlinear iterations = .*/set Max nonlinear iterations = $I /g" \
                   -e "s/set Linear solver tolerance =.*/set Linear solver tolerance = $i_LT/g" \
                   -e "s/set Nonlinear solver tolerance.*/set Nonlinear solver tolerance = $i_NLT/g" \
                   -e "s/set Linear solver A block tolerance.*/set Linear solver A block tolerance = $i_ABT/g" \
                   -e "s/    set Reference compressibility .*/    set Reference compressibility = $COMP/g" \
                   -e "s/set Max pre-Newton nonlinear iterations = .*/    set Max pre-Newton nonlinear iterations = $i_P/g" \
                   -e "s/set Use Newton failsafe = .*/set Use Newton failsafe = $i_UFS/g" \
                   -e "s/set Stabilization preconditioner = .*/set Stabilization preconditioner = $i_NSP/g" \
                   -e "s/set Stabilization velocity block = .*/set Stabilization velocity block = $i_NSA/g" \
                   -e "s/set SPD safety factor = .*/set SPD safety factor = $i_SF/g" \
                   -e "s/set Use deviator of strain-rate = .*/set Use deviator of strain-rate = $i_UDS/g" \
                   -e "s/set Max Newton line search iterations = .*/    set Max Newton line search iterations = $i_LS/g" \
                   -e "s/set Maximum linear Stokes solver tolerance =.*/set Maximum linear Stokes solver tolerance = $i_OS/g" \
                   -e "s/set Use Newton residual scaling method .*/    set Use Newton residual scaling method = $i_RSM/g" \
                  input.prm > "$infilename"
                  
                  nohup mpirun -np $processes $build_dir./aspect $infilename > $outfilename 2>$errorfilename
                  
                  grep "Relative nonlinear residual" $outfilename > $outplotfilename
                  
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
   done
  done
 done
done
