#!/bin/bash

# Note that it may take a long time to run all of these models, 
# depending on the resolution set in the input file. 
# I used 64 cores on a work station to generate the results, 
# which took about a day for refinement level 7.
# The models with layered convection are usually fast to run, 
# but each model with one large convection cell took about an
# hour.  

for rayleigh in "10000" "100000" "400000" "2000000"; do 

  for P in "0.4" "0.2" "0.1" "0.0" "-0.1" "-0.2" "-0.3" "-0.4" "-0.5" "-0.6" "-0.7" "-0.8"; do 

    rm current.prm

    k=$(echo "scale=5; 24603750/$rayleigh" | bc)
    gamma=$(echo "scale=5; $P/2.0 * 13500000" | bc)

    echo "subsection Material model" >> current.prm
    echo "  subsection Latent heat" >> current.prm
    echo "    set Thermal conductivity = $k" >> current.prm
    echo "    set Phase transition Clapeyron slopes = $gamma" >> current.prm
    echo "  end" >> current.prm
    echo "end" >> current.prm

    # write out statistics
    echo "set Output directory = Ra${rayleigh}_P${P}" >> current.prm
    echo "Starting Ra${rayleigh}_P${P}"
    cat christensen_yuen_phase_function.prm current.prm | mpirun -np 8 ../../build/aspect -- >log 2>log

  done
done
