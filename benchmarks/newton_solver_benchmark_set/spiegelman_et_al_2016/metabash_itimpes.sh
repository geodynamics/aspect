#!/bin/bash
version="c3"
declare -a grid=(4)
declare -a NSP=("none") # "symmetric" "none") #Use Newton stabilization preconditioner
declare -a NSA=("none") # "symmetric" "none") #Use Newton stabilization A block
declare -a UFS=("truefalse") #"false") # use failsafe
declare -a RSM=("truefalse") # "true") # use residual scaling
declare -a LS=(0-5-10) #(10) # line search
declare -a COMP=(0) # copmressibiltiy
declare -a phi=(30)
declare -a ST=("1e-20") # "1e-4")
declare -a UDS=("false" "true") #"true" #"false"
declare -a SF=("9e-1") # "9.9e-1") # "9.99e-1" "1")
declare -a AEW=("false") #"true" #"false"
processes=4

for i_grid in "${grid[@]}"
do
 for i_NSP in "${NSP[@]}"
 do
  for i_NSA in "${NSA[@]}"
  do
   for i_LS in "${LS[@]}"
   do
    for i_RSM in "${RSM[@]}"
    do
     for i_UFS in "${UFS[@]}"
     do
      for i_COMP in "${COMP[@]}"
      do
       for i_phi in "${phi[@]}"
       do
        for i_ST in "${ST[@]}"
        do
         for i_UDS in "${UDS[@]}"
         do
          for i_SF in "${SF[@]}"
          do
           for i_AEW in "${AEW[@]}"
           do
           output_name="$version""_g""$i_grid""_UFS""$i_UFS""_NSP-""$i_NSP""_NSA-""$i_NSA""_LS""$i_LS""_RSM""$i_RSM""_phi""$i_phi""_COMP""$i_COMP""_ST""$i_ST""_UDS""$i_UDS""_SF""$i_SF""_AEW""$i_AEW"
           log_file="auto_logs/""$output_name"".log"
           bash_auto_file="bash_autos/bash_""$output_name"".sh"
           echo "$output_name"
           sed \
            -e "s/declare -a grid.*/declare -a grid=($i_grid)/g" \
            -e "s/declare -a UFS.*/declare -a UFS=(\"false\" \"true\")/g" \
            -e "s/declare -a NSP.*/declare -a NSP=($i_NSP)/g" \
            -e "s/declare -a NSA.*/declare -a NSA=($i_NSA)/g" \
            -e "s/declare -a ST.*/declare -a ST=($i_ST)/g" \
            -e "s/declare -a UDS.*/declare -a UDS=($i_UDS)/g" \
            -e "s/declare -a SF.*/declare -a SF=($i_SF)/g" \
            -e "s/declare -a AEW.*/declare -a AEW=($i_AEW)/g" \
            -e "s/declare -a LS.*/declare -a LS=(0 5 10)/g" \
            -e "s/version=.*/declare -a version=\"$version\"/g" \
            -e "s/declare -a RSM.*/declare -a RSM=(\"false\" \"true\")/g" \
            -e "s/phi=.*/phi=\"$i_phi\"/g" \
            -e "s/COMP=.*/COMP=\"$i_COMP\"/g" \
            -e "s/processes=.*/processes=\"$processes\"/g" \
            bash.sh > $bash_auto_file
           
           nohup bash ./$bash_auto_file > $log_file 2>&1 &
           sleep 2 
           
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
