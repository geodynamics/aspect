# The script imports necessary libraries and sets up parameters for our temperature model.

import matplotlib.pyplot as plt
import numpy as np
import math

# Parameters
r_inner = 3481000  # Internal radius in meters
r_outer = 6371000  # External radius in meters
p_large = 400.0  # Amplitude of the large perturbations
p_small = 30.0  # Amplitude of the small perturbations
frequency_large = 3  # Frequency of the large perturbations
frequency_small = 30  # Frequency of the small perturbations
spread_angle_large_degrees = 90  # Spread angle of the large perturbations in degrees
start_angle_degrees = 330  # Starting angle for the perturbation in degrees

# The temperature_function_degrees_array function computes temperature values over a mesh grid of radii and angles, applying the perturbations and returning an array of temperature values.

def temperature_function_degrees_array(r, theta_degrees):
    # Adjust theta to shift the starting angle of the perturbation
    adjusted_theta_degrees = theta_degrees + start_angle_degrees
    adjusted_theta_degrees = adjusted_theta_degrees % 360

    # Initialize the temperature array
    temp = np.zeros_like(r)

    # Determine the large perturbation using a step function approach

    large_perturbation_mask = (adjusted_theta_degrees < spread_angle_large_degrees) | (adjusted_theta_degrees > 360 - spread_angle_large_degrees)
    temp[large_perturbation_mask] = p_large * np.sin(np.pi * (r[large_perturbation_mask] - r_inner) / (r_outer - r_inner)) * np.sin(frequency_large * np.radians(theta_degrees[large_perturbation_mask]))

    # Apply small perturbations where large perturbations are not applied
    small_perturbation_mask = ~large_perturbation_mask
    temp[small_perturbation_mask] = p_small * np.sin(np.pi * (r[small_perturbation_mask] - r_inner) / (r_outer - r_inner)) * np.sin(frequency_small * np.radians(theta_degrees[small_perturbation_mask]))

    return temp

# Here we generate data for visualization
theta = np.linspace(0, 360, 360)  # Angles from 0 to 360 degrees
radius = np.linspace(r_inner, r_outer, 100)  # Radii from r_inner to r_outer

# Here we create a grid for plotting

R, Theta = np.meshgrid(radius, theta)
Temp = temperature_function_degrees_array(R, Theta)

# After importing the libraries and defining our function, the script will plot how our temperature profile will look. The next step is to write the temperature perturbation in an ASCII format that ASPECT can read.
# The script includes a function, write_ascii_data, which is responsible for writing temperature data into an ASCII file that can be read by the ASPECT. The function takes arrays representing radial distances (radius), angular degrees (theta_degrees), temperature values (temperatures), and a filename (filename) where the data will be saved.
# The function opens a file with the given filename in write mode. All the data that follows will be written to this file.
# The line with the keyword "POINTS:" in the header of the ASCII file is crucial because it defines the dimensions of the data grid for the ASPECT simulation program. The numbers following "POINTS:" represent the size of the grid in each dimension:
# 100 represents the number of points along the radial dimension (nx). In this case, it means there are 100 radial divisions or layers from the inner radius r_inner to the outer radius r_outer.
# 360 represents the number of points along the angular dimension (ny). This indicates that the angular data is divided into 360 segments, corresponding to the degrees in a full circle.
# 1 represents the number of points along the third dimension (nz). In this context, it implies that the data is essentially two-dimensional, with no depth or third dimension to vary over. This would be applicable for a surface or a shell where only the radius and the angle change.
# The structure indicated by "POINTS: 100 360 1" tells ASPECT how to read the subsequent data lines. Each line of data will correspond to a point in this grid, and ASPECT will expect to read 100 x 360 x 1 lines of data following the headers. The arrangement of data points in the file must match this grid layout for ASPECT to interpret the file correctly.
#In the write_ascii_data function, the nested for loops iterate over the angular positions (theta_degrees) and radial positions (radius), writing out the temperature (temp) at each point. The formatting ensures that each value occupies a specific width in the file, aligning the data into neat columns. This formatting is not just for aesthetics; it is essential for ensuring that ASPECT can correctly parse each value, especially when the data is read by automated scripts or tools that expect a consistent column width. Without such structured formatting, there could be parsing errors, leading to incorrect simulation results or program failures.
# The function iterates through the provided arrays of angles and radii to write the temperature data for each grid point.
# For each point, it formats the radius (r), angle (phi), and temperature (temp) with improved formatting to ensure numerical precision and readability. The formatting specifies that radius should have one decimal place, angle should have four decimal places (or be an integer if it is a whole number), and temperature should have one decimal place (or be an integer if it is a whole number).
# The data is written with aligned columns, ensuring that each row lines up correctly for the radius, angle, and temperature values. This is crucial for ASPECT to correctly parse the file. Each value is given a specific column width with left alignment.
# At the end, we convert degrees to radians, since ASPECT only works with radians. Then it creates the file perturbation.txt



def write_ascii_data(radius, theta_degrees, temperatures, filename):
    with open(filename, 'w') as file:
        file.write("# Test data for ascii data initial conditions.\n")
        file.write("# Only next line is parsed in format: [nx] [ny] [nz] because of keyword \"POINTS:\"\n")
        file.write("# POINTS: 100 360 1\n")
        file.write("# Columns: r phi      temperature [K]\n")

        for j in range(len(theta_degrees)):
            for i in range(len(radius)):
                r = radius[i]
                phi = theta_degrees[j]
                temp = temperatures[j, i]

                # Improved formatting
                r_str = f"{r:.1f}"
                phi_str = f"{phi:.4f}" if not phi.is_integer() else f"{int(phi)}"
                temp_str = f"{temp:.1f}" if not temp.is_integer() else f"{int(temp)}"

                # Column alignment with offset for the second column
                file.write(f"{r_str:<10} {phi_str:<11}{temp_str}\n")


# Conversion of polar coordinates to Cartesian coordinates for plotting
X = R * np.cos(np.radians(Theta))
Y = R * np.sin(np.radians(Theta))


# Plotting
plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, Temp, levels=50, cmap='viridis')
plt.colorbar(contour, label='Temperature')
plt.title('Temperature distribution')
plt.xlabel('Coordinate X (m)')
plt.ylabel('Coordinate Y (m)')
plt.axis('equal')
plt.show()


# Use the function write_ascii_data to save data
write_ascii_data(radius, theta, Temp, 'perturbation_ascii.txt')

# In this step, we will transform the angle measurements from degrees into radians. It's crucial to make this conversion because ASPECT is compatible with ASCII files that contain radian values only, not degrees.

import math

def convert_to_radians(degrees):
    """Convert degrees to radians."""
    return degrees * (math.pi / 180)

def convert_file_angles(input_file_path, output_file_path):
    """Convert the angles in the second column of the file from degrees to radians."""
    with open(input_file_path, 'r') as file:
        content = file.readlines()

    converted_lines = content[:4]  # Keep the header lines unchanged
    for line in content[4:]:
        columns = line.split()
        if len(columns) == 3:
            radius, angle_deg, temperature = columns
            angle_rad = convert_to_radians(float(angle_deg))
            converted_line = f"{radius}  {angle_rad:.8f}  {temperature}\n"
            converted_lines.append(converted_line)
        else:
            converted_lines.append(line)

    with open(output_file_path, 'w') as file:
        file.writelines(converted_lines)

# Example usage
input_file_path = 'perturbation_ascii.txt'  # Replace with your input file path
output_file_path = 'perturbation_ascii.txt'  # Replace with your desired output file path

convert_file_angles(input_file_path, output_file_path)
