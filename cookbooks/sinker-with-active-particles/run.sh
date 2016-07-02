# aspect executable 
ASPECT="/Users/jbnaliboff/Software/aspect/aspect/./aspect"

# Model name
name="sinker-with-active-particles"

mpirun -np 1 $ASPECT $name".prm" &> "stdout_"$name &
