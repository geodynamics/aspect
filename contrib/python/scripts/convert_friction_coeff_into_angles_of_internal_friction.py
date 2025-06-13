# This script can be used to convert a given list of friction coefficients
# into angles of internal friction. This is necessary, because 
# most ASPECT material models expect angles of internal friction as input.
#
# The friction coefficient is defined as tan(angle_of_internal_friction)
# where angle_of_internal_friction is in radians.
# To convert the friction coefficient to angle of internal friction in degrees:
# angle_of_internal_friction = arctan(friction_coefficient) * 180 / pi

import numpy as np

friction_coefficient_list = [0.6,0.8,1.3] # List of example friction coefficients. Change these to reflect the friction coefficients appropriate for your composition.

for i in range (len(friction_coefficient_list)):
    angle_of_internal_friction = np.arctan(friction_coefficient_list[i]) * 180 / np.pi
    print("Friction coefficient: " + str(friction_coefficient_list[i]) + ", Angle of internal friction: " + str(angle_of_internal_friction))

    
