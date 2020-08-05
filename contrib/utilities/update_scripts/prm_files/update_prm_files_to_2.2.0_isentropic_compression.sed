# A script for the stream editor sed to update .prm files for ASPECT
# to account for the new name of the isentropic compression formulation.
# This formulation was previously named isothermal compression, which
# was misleading and wrong. 

s/isothermal compression/isentropic compression/g
