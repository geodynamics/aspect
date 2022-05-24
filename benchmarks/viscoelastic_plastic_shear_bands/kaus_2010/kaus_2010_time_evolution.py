# Based on the following script:
# benchmarks/viscoelastic_plastic_shear_bands/gerya_2019/gerya_2019_analysis.py
# Help: python3 kaus_2010_time_evolution.py --help
import argparse
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import vtk as vtk
from vtk.util import numpy_support


def main():
  parser = argparse.ArgumentParser(
    description='Create a strain rate profile from the output directory of the aspect',
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

  plot_timesteps(
    solution_folder=output_dir / 'solution',
    plot_file=plot_file,
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


def plot_timesteps(solution_folder, plot_file, interesting_times, xsize=8, ysize=2, use_log=True, vmin=None, vmax=None):
  fig = plt.figure(
    figsize=(xsize, ysize*len(interesting_times)), # in inch
    # facecolor='white',
  )
  fig.suptitle(
    'Strain Rate Second Invariant in $\log_{10}$ scale (1/s)',
    x=0.5, y=0.98,
    fontsize=14
  )
  fig.supxlabel('Horizontal Position (m)')
  fig.supylabel('Vertical Position (m)')

  gs = fig.add_gridspec(
    nrows=len(interesting_times), ncols=1,
    left=0.125, right=0.9,
    bottom=0.11, top=0.88,
    hspace=0.2, wspace=0.2,
  )

  for i, time in enumerate(interesting_times):
    # -------------------------------------------------------------------------
    # Load the data
    # -------------------------------------------------------------------------
    # Find the last pvtu file in the solution folder
    pvtu_file = sorted(solution_folder.glob(f"solution-{time:05d}.pvtu"))[0]

    x, y, field = get_data(pvtu_file, 'strain_rate')


    # -------------------------------------------------------------------------
    # Plot: Strain Rate Field
    # -------------------------------------------------------------------------
    ax = fig.add_subplot(
      gs[i, 0],
      aspect='equal',
      # facecolor='white',
      title=f"Timestep: {time}",
      # xlabel='Horizontal Position (m)',
      xlim=(np.min(x), np.max(x)),
      xticks=[] if i < len(interesting_times) - 1 else range(0,45000, 5000), # Hide if empty
      # xticklabels=[], # Hide if empty
      # ylabel='Vertical Position (m)',
      ylim=(np.min(y), np.max(y)),
      # yticks=[], # Hide if empty
      # yticklabels=[], # Hide if empty
    )
    sc = ax.scatter(
      x, y,
      s=4, marker='s',
      c=np.log10(field), edgecolors='none',
      vmin=np.log10(vmin) if vmin else None,
      vmax=np.log10(vmax) if vmax else None,
      rasterized=True,
    )
    if not (vmin and vmax):
      fig.colorbar(
        sc, ax=ax,
        # label='',
        # ticks=[], # Hide if empty
        location='right',
        orientation='vertical',
        fraction=0.1,
        aspect=30,
      )

  if vmin and vmax:
    fig.colorbar(
      sc, ax=fig.get_axes(),
      # label='',
      # ticks=[], # Hide if empty
      location='right',
      orientation='vertical',
      fraction=0.1,
      aspect=30,
    )

  plot_file.parent.mkdir(parents=True, exist_ok=True)
  plt.savefig(plot_file)


if __name__ == "__main__":
  main()
