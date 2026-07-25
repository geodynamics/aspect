import numpy as np
import matplotlib.pyplot as plt
import os
from vtk import vtkXMLUnstructuredGridReader
from vtk.util.numpy_support import vtk_to_numpy


def extract_vtu_data(soln_dir, time_step, variable_names):
    '''
    A function which uses the vtk python package to extract the quadrature points and associated scalar field values from the .vtu files output by ASPECT.

    soln_dir:       The absolute path to the solution directory output by ASPECT, which contains all the .vtu files.
    time_step:      The number on the solution file output to be loaded. ASPECT outputs .vtu files in the form of solution-XXXXX-YYYY.vtu,
                    where X corresponds to the time_step, and Y corresponds to the process number.
    variable_names: A list of variable names which correspond to the field names to be extracted from the ASPECT .vtu files.
    '''
    # Sort the files in the solution directory, and create empty arrays for storing the quadrature points
    # and create a list of empty lists equal to the number of desired variables to store the scalar fields.
    solutions   = np.sort(os.listdir( soln_dir ))
    coordinates = []
    extracted_variables = [[] for _ in range(len(variable_names))]

    # Iterate over all files in the solutions directory, check to make sure the files are .vtu files, and
    # that they are the appropriate time step.
    for file in solutions:
        if file.endswith('.vtu') and file[9:14] == str(time_step).zfill(5):
            # Absolute path to the current .vtu file
            file_path = os.path.join(soln_dir, file)

            # These next lines create an unstructured grid reader, which parses the .vtu file
            reader = vtkXMLUnstructuredGridReader()
            reader.SetFileName(file_path)
            reader.Update()
            data = reader.GetOutput()
            points = data.GetPoints()

            # Extract the quadrature points, and iterate over the desired scalar fields, adding them
            # to the storage array.
            x = vtk_to_numpy(points.GetData())
            coordinates.append(x)

            for k in range(len(extracted_variables)):
                if variable_names[k] == "stress_xx":
                    stresses = vtk_to_numpy(data.GetPointData().GetArray("stress"))
                    extracted_variables[k].extend(stresses[:, 0])

                elif variable_names[k] == "stress_yy":
                    stresses = vtk_to_numpy(data.GetPointData().GetArray("stress"))
                    extracted_variables[k].extend(stresses[:, 4])

                elif variable_names[k] == "shear_stress":
                    shears = vtk_to_numpy(data.GetPointData().GetArray("shear_stress"))
                    extracted_variables[k].extend(shears[:, 1])

                else:
                    # Extract the analytical solution from the .vtu files
                    extracted_variables[k].extend(vtk_to_numpy(data.GetPointData().GetArray(variable_names[k])))

    # Concatenate the points, creating a numpy array in the process
    coordinates = np.concatenate(coordinates, axis=0)

    # Output the extracted scalar fields as a numpy array.
    return coordinates, np.array(extracted_variables).T

def extract_data_along_line(coordinates, extracted_variables):
    '''
    A function which extracts the analytical solution along a line in the domain, and outputs it as a numpy array.
    '''
    # Create an array of x values along the line, and an empty array to store the analytical solution.
    vertical_profiles = np.array([300e3])
    surface_profiles  = np.array([np.max(coordinates[:, 1])])

    surface_indices   = np.where(coordinates[:, 1] == surface_profiles[0])[0]
    vertical_indices  = np.where(coordinates[:, 0] == vertical_profiles[0])[0]

    surface_x         = coordinates[surface_indices, 0]
    surface_stresses  = extracted_variables[surface_indices, :]

    vertical_y        = coordinates[vertical_indices, 1]
    vertical_stresses = extracted_variables[vertical_indices, :]

    return surface_x, surface_stresses, vertical_y, vertical_stresses

def plot_comparison(coordinates, extracted_variables, variable_names):
    '''
    A function which plots the analytical solution and the ASPECT solution for each of the desired variables.

    coordinates:        The quadrature points extracted from the .vtu files.
    extracted_variables: The scalar fields extracted from the .vtu files.
    variable_names:     A list of variable names which correspond to the field names to be extracted from the ASPECT .vtu files.
    '''
    # Iterate over all desired variables, and plot them against the analytical solution.
    surface_x, surface_stresses, vertical_y, vertical_stresses = extract_data_along_line(coordinates, extracted_variables)
    depth = 300e3 - vertical_y

    analytic_color = "mediumblue"
    aspect_color   = "orange"


    axis_label_fontsize = 14
    title_fontsize      = 20

    # Plot 2 rows of 3 columns. The first row is the surface stresses, and the second row is the vertical stresses.
    plt.figure(dpi=160, figsize=(16, 8))
    plt.subplot(2, 3, 1)

    plt.plot(surface_x / 1e3, surface_stresses[:, 3] / 1e6, label='ASPECT', lw=5, color=aspect_color)
    plt.plot(surface_x / 1e3, surface_stresses[:, 0] / 1e6, label='Analytic', lw=3, color=analytic_color)

    plt.xlabel('x - km', fontsize=axis_label_fontsize)
    plt.ylabel("$\\sigma_{xx}$ - MPa", fontsize=axis_label_fontsize)
    plt.legend()


    plt.subplot(2, 3, 2)
    plt.title("Stresses at the Surface", fontsize=title_fontsize)
    plt.plot(surface_x / 1e3, surface_stresses[:, 4] / 1e6, label='ASPECT', lw=5, color=aspect_color)
    plt.plot(surface_x / 1e3, surface_stresses[:, 1] / 1e6, label='Analytic', lw=3, color=analytic_color)

    plt.xlabel('x - km', fontsize=axis_label_fontsize)
    plt.ylabel("$\\sigma_{xy}$ - MPa", fontsize=axis_label_fontsize)

    plt.subplot(2, 3, 3)
    plt.plot(surface_x / 1e3, surface_stresses[:, 5] / 1e6, label='ASPECT', lw=5, color=aspect_color)
    plt.plot(surface_x / 1e3, surface_stresses[:, 2] / 1e6, label='Analytic', lw=3, color=analytic_color)

    plt.xlabel('x - km', fontsize=axis_label_fontsize)
    plt.ylabel("$\\sigma_{yy}$ - MPa", fontsize=axis_label_fontsize)

    plt.subplot(2, 3, 4)
    plt.plot(vertical_stresses[:, 3] / 1e6, depth / 1e3, label='ASPECT', lw=5, color=aspect_color)
    plt.plot(vertical_stresses[:, 0] / 1e6, depth / 1e3, label='Analytic', lw=3, color=analytic_color)
    plt.gca().invert_yaxis()

    plt.xlabel("$\\sigma_{xx}$ - MPa", fontsize=axis_label_fontsize)
    plt.ylabel('depth - km', fontsize=axis_label_fontsize)

    plt.subplot(2, 3, 5)
    plt.title("Stresses with Depth Beneath the Load", fontsize=title_fontsize)
    plt.plot(vertical_stresses[:, 4] / 1e6, depth / 1e3, label='ASPECT', lw=5, color=aspect_color)
    plt.plot(vertical_stresses[:, 1] / 1e6, depth / 1e3, label='Analytic', lw=3, color=analytic_color)
    plt.gca().invert_yaxis()

    plt.xlabel("$\\sigma_{xy}$ - MPa", fontsize=axis_label_fontsize)
    plt.ylabel('depth - km', fontsize=axis_label_fontsize)

    plt.subplot(2, 3, 6)
    plt.plot(vertical_stresses[:, 5] / 1e6, depth / 1e3, label='ASPECT', lw=5, color=aspect_color)
    plt.plot(vertical_stresses[:, 2] / 1e6, depth / 1e3, label='Analytic', lw=3, color=analytic_color)
    plt.gca().invert_yaxis()

    plt.xlabel("$\\sigma_{yy}$ - MPa", fontsize=axis_label_fontsize)
    plt.ylabel('depth - km', fontsize=axis_label_fontsize)

    plt.tight_layout()
    plt.savefig(f'surface_and_depth_stresses.png', dpi=300, bbox_inches='tight')
    plt.close()


variable_names = ["flamant_sigma_xx", "flamant_sigma_xy", "flamant_sigma_yy", "stress_xx", "shear_stress", "stress_yy"]
soln_dir       = "../elastic_line_load/solution/"
time_step      = 1

coordinates, extracted_variables = extract_vtu_data(soln_dir, time_step, variable_names)
plot_comparison(coordinates, extracted_variables, variable_names)
