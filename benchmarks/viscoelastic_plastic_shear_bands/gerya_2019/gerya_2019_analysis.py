# Python script to load the solution data into a
# numpy array, interpolate the data to a new
# uniform grid, and plot the results.
# From a terminal, execute the script with
#   python gerya_2019_analysis.py
# Within a python interpreter (e.g., a Jupyter notebook), 
# the script can be executed with the command:
#   exec(open("gerya_2019_analysis.py").read())
# This script was tested with the following 
# package versions:
#   python     - 3.7.7
#   numpy      - 1.19.1
#   scipy      - 1.5.2
#   vtk        - 8.2.0
#   matplotlib - 3.3.1

# Load modules
import numpy as np
import matplotlib.pyplot as plt
import vtk as vtk; from vtk.util import numpy_support
from scipy.interpolate import griddata

# Relative (from current directory) path to the
# solution data folder.
vtk_directory = 'output_gerya_2019_vp/solution/'

# Time step to analyze. This number corresponds to 
# the numbering of .pvtu files. In this instance
# '00000.0213' is the last nonlinear iteration
# of the first time step, but note that the 
# number of nonlinear iterations may vary slightly
# between machines, or significantly when additional
# parameters are altered.
pvtu_number = '00000.0213'

# Load vtu data (pvtu directs to vtu files)
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName(vtk_directory + 'solution-' + pvtu_number + '.pvtu')
reader.Update()

# Get the coordinates of nodes in the mesh
nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

# Convert nodal vtk data to a numpy array
nodes_numpy_array = vtk.util.numpy_support.vtk_to_numpy(nodes_vtk_array)

# Extract x, y and z coordinates from numpy array 
x,y = nodes_numpy_array[:,0] , nodes_numpy_array[:,1]

# Determine the number of scalar fields contained in the .pvtu file
number_of_fields = reader.GetOutput().GetPointData().GetNumberOfArrays()

# Determine the name of each field and place it in an array.
field_names = []
for i in range(number_of_fields):
  field_names.append(reader.GetOutput().GetPointData().GetArrayName(i))

# Determine the index of the field strain_rate
idx = field_names.index("strain_rate")

# Extract values of strain_rate
field_vtk_array = reader.GetOutput().GetPointData().GetArray(idx)
strain_rate     = numpy_support.vtk_to_numpy(field_vtk_array)

# Generate coordinate arrays for data interpolation based on model dimensions
# and grid spacing. The coordinate arrays will have the double the resolution
# as the grid spacing.
res = x[1] - x[0]
xi = np.linspace(np.min(x), np.max(x), int((np.max(x) - np.min(x))/res)*2 + 1 )
yi = np.linspace(np.min(y), np.max(y), int((np.max(y) - np.min(y))/res)*2 + 1 )

# Generate data interpolation grid. 
X, Y = np.meshgrid(xi,yi)

# Interpolate strain rate onto new grid
R = griddata((x,y), strain_rate, (X,Y), method='cubic')

# Height of strain-rate profile (model height * 5/16)
profile_height = np.max(y) * 5./16. 

# Find all indices corresponding the profile
inds = np.where(Y[:]==profile_height)

# New array for x and strain rate along profile
x_r = np.concatenate((X[inds].reshape(-1,1),R[inds].reshape(-1,1)),axis=1)

# Plot strain-rate field
plt.figure()
fig, (ax0) = plt.subplots(1,1)
c = ax0.pcolor(X,Y,R,shading='auto')
ax0.set_aspect('equal', 'box')
ax0.set_xlabel('Horizontal Position (m)')
ax0.set_ylabel('Vertical Position (m)')
ax0.set_title('Strain Rate Second Invariant (1/s)', fontsize=10)
fig.colorbar(c, ax=ax0)
plt.savefig("strain_rate_field.png", dpi=300)
plt.close()

# Plot strain rate profile
plt.figure()
fig, (ax0) = plt.subplots(1,1)
c = ax0.plot(x_r[:,0],np.log10(x_r[:,1]))
ax0.set_xlabel('Horizontal Position (m)')
ax0.set_ylabel('Log10 of Strain Rate Second Invariant (1/s)')
ax0.set_title('Strain Rate Profile at y = Model_Height * 5/16', fontsize=10)
plt.savefig("strain_rate_profile.png", dpi=300)
plt.close()
