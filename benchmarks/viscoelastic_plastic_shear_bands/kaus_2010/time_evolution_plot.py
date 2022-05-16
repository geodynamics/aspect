# Based on the following script:
# benchmarks/viscoelastic_plastic_shear_bands/gerya_2019/gerya_2019_analysis.py
# Usage: python3 time_evolution_plot.py -i <OUTPUT_FOLDER_FROM_THE_PRM> -o <PATH_TO_OUTPUT_DIR_FOR_PLOTS>
import argparse
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import vtk as vtk
from vtk.util import numpy_support


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "-i", "--input",
    type=str,
    required=True,
  )
  parser.add_argument(
    "-o", "--output",
    type=str,
    required=True,
  )
  args = parser.parse_args()

  input_dir = Path(args.input)
  output_dir = Path(args.output)

  plot_timesteps(
    solution_folder=input_dir / 'solution',
    output_dir=output_dir,
    interesting_times=[6, 12, 18],
    xsize=8,
    ysize=2,
    use_log=True,
    vmin=1e-15,
    vmax=1e-13,
  )


def get_data(file_path, field_name):
  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(file_path)
  reader.Update()

  # Get the coordinates of nodes in the mesh
  nodes_vtk_array = reader.GetOutput().GetPoints().GetData()

  # Convert nodal vtk data to a numpy array
  nodes_numpy_array = vtk.util.numpy_support.vtk_to_numpy(nodes_vtk_array)

  # Extract x, y and z coordinates from numpy array 
  x, y = nodes_numpy_array[:,0] , nodes_numpy_array[:,1]

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


def reshape_data(x, y, field_array, res=None):
  if (res is None):
    res = x[1] - x[0]
    xi = np.linspace(np.min(x), np.max(x), int((np.max(x) - np.min(x))/res)*2 + 1 )
    yi = np.linspace(np.min(y), np.max(y), int((np.max(y) - np.min(y))/res)*2 + 1 )
  else:
    resolution_x = res
    resolution_y = int(resolution_x * (np.max(y) - np.min(y)) / (np.max(x) - np.min(x)))
    xi = np.linspace(np.min(x), np.max(x), resolution_x)
    yi = np.linspace(np.min(y), np.max(y), resolution_y)

  # Generate data interpolation grid
  X, Y = np.meshgrid(xi, yi)

  # Interpolate strain rate onto new grid
  R = griddata((x,y), field_array, (X,Y), method='linear')

  return X, Y, R


def generate_field_map(x, y, field_array):
  fig = plt.figure(figsize=(8, 2))

  ax = fig.add_subplot(1, 1, 1)

  c = ax.pcolor(x, y, field_array, shading='auto')
  ax.set_aspect('equal', 'box')
  ax.set_xlabel('Horizontal Position (m)')
  ax.set_ylabel('Vertical Position (m)')
  ax.set_title('Strain Rate Second Invariant (1/s)', fontsize=10)
  fig.colorbar(c, ax=ax)

  plt.show()


def plot_timesteps(solution_folder, output_dir, interesting_times, xsize=8, ysize=2, use_log=True, vmin=None, vmax=None):
  fig, axes = plt.subplots(len(interesting_times), 1, sharex=True)
  fig.set_size_inches((xsize, ysize*len(interesting_times)))

  for i, time in enumerate(interesting_times):
    pvtu_file = sorted(solution_folder.glob(f"solution-{time:05d}.pvtu"))[0]

    x, y, field = get_data(pvtu_file, 'strain_rate')
    x, y, field = reshape_data(x, y, field, res=100)

    ax = axes[i]

    if use_log:
      c = ax.pcolor(x, y, np.log10(field), shading='auto', vmin=np.log10(vmin) if vmin else None, vmax=np.log10(vmax) if vmax else None)
    else:
      c = ax.pcolor(x, y, field, shading='auto', vmin=vmin, vmax=vmax)

    ax.set_aspect('equal', 'box')
    ax.label_outer()

    ax.set_title(f"Timestep: {time}")

    if not (vmin and vmax):
      fig.colorbar(c, ax=ax)

  axes[1].set_ylabel('Vertical Position (m)')
  axes[2].set_xlabel('Horizontal Position (m)')

  fig.tight_layout()

  fig.suptitle('Strain Rate Second Invariant log10(1/s)', x=0.5, y=1.025)

  if vmin and vmax:
    fig.colorbar(c, ax=axes.ravel().tolist())

  output_dir.mkdir(parents=True, exist_ok=True)
  plt.savefig(output_dir / 'strain_rate_fields_evolution.pdf',  format='pdf', bbox_inches='tight')


if __name__ == "__main__":
  main()
