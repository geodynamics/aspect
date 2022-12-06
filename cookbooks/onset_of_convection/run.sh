#!/bin/bash

# Number of processors to run on
NP=1
filename="onset-convection-data.csv"
mantle_thickness=3000000
pi=3.1415926536

rm $filename
echo "# Viscosity (Pa s) DeltaT (K) velocity (m/yr)" >> $filename

echo "----Running models now----"
for viscosity in "1e24" "3.3e24" "1e25" "3.3e25" "1e26" "3.3e26" "1e27"
do
  for DeltaT in "10" "30" "100" "300" "1000" "3000" "10000"
  do  
    echo "Running model with a viscosity of $viscosity Pa s and a temperature jump of $DeltaT:"
    cp onset_of_convection.prm temp.prm
    name=$viscosity"_"$DeltaT
    echo "set Output directory = output_convection/convection_$name" >> temp.prm
    echo "subsection Termination criteria" >> temp.prm
    echo "  set End step                             = 2" >> temp.prm
    echo "end" >> temp.prm

    model_width=$(expr $mantle_thickness*3.1415926536 | bc)
    echo "subsection Geometry model" >> temp.prm
    echo "  subsection Box" >> temp.prm
    echo "    set Y extent = $mantle_thickness" >> temp.prm
    echo "    set X extent = $model_width" >> temp.prm
    echo "  end" >> temp.prm
    echo "end" >> temp.prm

    echo "subsection Initial temperature model" >> temp.prm
    echo "  subsection Function" >> temp.prm
    echo "    set Function constants  = p=1, L=$model_width, H=$mantle_thickness, pi=3.1415926536, k=2" >> temp.prm
    echo "    set Function expression = $DeltaT * (1.0-z/H) - p*cos(k*pi*x/L)*sin(pi*z/H)" >> temp.prm
    echo "  end" >> temp.prm
    echo "end" >> temp.prm

    echo "subsection Boundary temperature model" >> temp.prm
    echo "  subsection Box" >> temp.prm
    echo "    set Bottom temperature = $DeltaT" >> temp.prm
    echo "  end" >> temp.prm
    echo "end" >> temp.prm

    echo "subsection Material model" >> temp.prm
    echo "  subsection Simple model" >> temp.prm
    echo "    set Viscosity                     = $viscosity" >> temp.prm
    echo "  end" >> temp.prm
    echo "end" >> temp.prm

    echo "subsection Postprocess" >> temp.prm
    echo "  subsection Visualization" >> temp.prm
    echo "    set List of output variables = material properties" >> temp.prm
    echo "    subsection Material properties" >> temp.prm
    echo "      set List of material properties = density, thermal expansivity, specific heat, viscosity, thermal conductivity" >> temp.prm
    echo "    end" >> temp.prm
    echo "  end" >> temp.prm
    echo "end" >> temp.prm

    printf "$viscosity $DeltaT " >> $filename
    mpirun -n $NP $1 temp.prm | tail -41 | gawk '/RMS, max velocity/ {printf "%s ", $4}' | sed 's/,//g' >> $filename
    printf "\n" >> $filename
    rm -f temp.prm
  done
  printf "\n" >> $filename
done
