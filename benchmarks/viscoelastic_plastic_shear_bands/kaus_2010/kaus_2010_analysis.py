# Based on the following script:
# benchmarks/viscoelastic_plastic_shear_bands/gerya_2019/gerya_2019_analysis.py
# Help: python3 kaus_2010_analysis.py --help
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import vtk as vtk
from vtk.util import numpy_support



def get_data(pvtu_file, field_name):
  # Originally created by John Naliboff
  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(pvtu_file)
  reader.Update()

  # Get the coordinates of nodes in the mesh
  nodes_vtk_array = reader.GetOutput().GetPoints().GetData()

  # Convert nodal vtk data to a numpy array
  nodes_numpy_array = vtk.util.numpy_support.vtk_to_numpy(nodes_vtk_array)

  # Extract x, y and z coordinates from numpy array 
  x, y = nodes_numpy_array[:,0], nodes_numpy_array[:,1]

  # Determine the number of scalar fields contained in the .pvtu file
  number_of_fields = reader.GetOutput().GetPointData().GetNumberOfArrays()

  # Determine the name of each field and place it in an array.
  field_names = []
  for i in range(number_of_fields):
    field_names.append(reader.GetOutput().GetPointData().GetArrayName(i))

  # Determine the index of the field strain_rate
  idx = field_names.index(field_name)

  # Extract values of strain_rate
  field_vtk_array = reader.GetOutput().GetPointData().GetArray(idx)
  field_array     = numpy_support.vtk_to_numpy(field_vtk_array)

  return x, y, field_array



def main():
  parser = argparse.ArgumentParser(
    description='Create a strain rate field and profile from the output directory of the aspect',
    epilog='Author: Marcel Saaro, John Naliboff'
  )
  parser.add_argument(
    "-o", "--output-dir",
    type=str,
    required=True,
    help='Path of the output dir in the prm file',
  )
  parser.add_argument(
    "-p", "--plot-file",
    type=str,
    required=True,
    help='Path of the plot file. None existing folders in the path will be created.',
  )
  args = parser.parse_args()

  output_dir = Path(args.output_dir)
  plot_file = Path(args.plot_file)

  # -------------------------------------------------------------------------
  # Load the data
  # -------------------------------------------------------------------------
  solution_folder = Path(output_dir) / 'solution'
  # Find the last pvtu file in the solution folder
  pvtu_file = sorted(solution_folder.glob('solution-*.pvtu'))[-1]

  x, y, strain_rate = get_data(pvtu_file, 'strain_rate')


  # -------------------------------------------------------------------------
  # Generate the Profile
  # -------------------------------------------------------------------------
  profile_indx = np.argwhere(abs(y  - 17./25.*np.max(y)) < 1e2)

  df = pd.DataFrame({
    'x' : x[profile_indx].ravel(),
    'profile' : strain_rate[profile_indx].ravel(),
  })
  df = df.groupby('x').mean()


  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  fig = plt.figure(
    figsize=(12, 4),
  )

  gs = fig.add_gridspec(
    nrows=1, ncols=1,
    left=0.125, right=0.9,
    bottom=0.11, top=0.88,
    hspace=0.2, wspace=0.2,
  )

  # Strain Rate Profile
  ax = fig.add_subplot(
    gs[0, 0],
    title='Strain Rate Profile at y = Model_Height * 5/16',
    xlabel='Horizontal Position (m)',
    ylabel='Log10 of Strain Rate Second Invariant (1/s)',
  )
  df['profile'].plot(ax=ax)

  # Save the figure
  plot_file.parent.mkdir(parents=True, exist_ok=True)
  plt.savefig(plot_file)



if __name__ == "__main__":
  main()
