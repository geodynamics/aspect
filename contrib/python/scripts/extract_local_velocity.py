# This script uses Paraview's pvpython to take output from a global mantle
# convection model and extract velocities along user-specified boundaries. 
# Next, the script uses the extracted velocities to create ASCII files which
# can be applied as boundary conditions in regional models to account for
# far-field effects outside of the regional model boundaries. This is achieved
# by specifying the bounds of the regional model, slicing through the 
# global model along planes which coincide with the location of the
# regional model boundaries, and saving the velocities at each point
# into an ASCII file. The regional model is expected to be a 3D spherical
# chunk, and the global model is expected to be a 3D spherical shell.

# The user must specify the following parameters in the script:
# 1. The location of the .pvd file for a global convection model
# 2. The output directory where the ASCII files will be saved
# 3. The refinement level of the global model 
# 4. The radial resolution of the regional model
# 5. The lateral resolution of the regional model
# 6. The maximum and minimum radius, latitude, and  longitude for the regional model.

# The default settings inside this script are for an example applying this script to
# the global model presented in the S2ORTS cookbook, and applying the extracted 
# velocities to a regional model that spans from 20 degrees latitude to 50 degrees
# latitude, 190 degrees longitude to 230 degrees longitude, and from a
# radius of 5070 km to a radius of 6370 km. Before running this script, 
# make sure that you run the S2ORTS cookbook and generate the solution.pvd file.

# This script defines 3 functions
# 1. spherical_to_cartesian, which converts from spherical coordinates to Cartesian
#    coordinates
# 2. cartesian_to_spherical, which converts from Cartesian coordinates to spherical
#    coordinates
# 3. slice_plane_calculator, which determines the normal to the plane that defines
#    the east and west model boundaries in the regional chunk.

# Import packages
import numpy as np
import pandas as pd
import os
from paraview import simple
from paraview import servermanager

def spherical_to_cartesian(radii, latitudes, longitudes, for_slice):
    """
    Converts from spherical coordinates to Cartesian coordinates.
    radii: the radius, in m
    latitudes: the latitude, in degrees ranging from -90 to 90
    longitudes: the longitude, in degrees ranging from 0 to 360
    for_slice: boolean, if false the provided points are directly 
    converted into Cartesian coordinates. If true, radii, latitudes,
    and longitudes must be arrays, and the function returns a uniform 
    structured grid defining the slice.
    Returns x, y, z, in m
    """

    if for_slice:
        cartesian_coordinates = []
        for lat in latitudes:
            for lon in longitudes:
                for r in radii:
                    x = r * np.sin(np.deg2rad(90 - lat)) * np.cos(np.deg2rad(lon))
                    y = r * np.sin(np.deg2rad(90 - lat)) * np.sin(np.deg2rad(lon))
                    z = r * np.cos(np.deg2rad(90 - lat))
    
                    if lon == 0:
                        y = 0
                    elif lon == np.pi/2:
                        x = 0
                    if lat == np.pi/2:
                        z = 0          
                    cartesian_coordinates.append([x, y, z])
    
        return np.array(cartesian_coordinates)
        
    else:
        x = radii * np.sin(np.deg2rad(90 - latitudes)) * np.cos(np.deg2rad(longitudes))
        y = radii * np.sin(np.deg2rad(90 - latitudes)) * np.sin(np.deg2rad(longitudes))
        z = radii * np.cos(np.deg2rad(90 - latitudes))
    
        return np.array([x, y, z])
    
def cartesian_to_spherical(x, y, z):
    """
    Takes an x, y, z point and converts it to spherical coordinates. 
    Returns r in m, latitude in degrees, and longitude, ranging from 0 to 360, in degrees 
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    latitude = 90 - np.rad2deg( np.arccos( z / (np.sqrt(x**2 + y**2 + z**2)) ) )
    longitude =  np.sign(y) * np.rad2deg(np.arccos( x / np.sqrt(x**2 + y**2) ))
    longitude[np.where(longitude < 0)] = longitude[np.where(longitude < 0)] + 360
    longitude[np.where(longitude == 0)] = 180
    
    return r, latitude, longitude

def slice_plane_calculator(boundary_name, radius_bounds, latitude_bounds, longitude_bounds):
    """
    Calculates the normal of a plane which is used for slicing the global models. This is
    achieved by defining three points on either the east or west model boundary using the 
    values provided by radius_bounds, latitude_bounds, and longitude_bounds.
    boundary_name: the name of the model boundary
    radius_bounds: the maximum and minimum radius of the regional models
    latitude_bounds: the maximum and minimum latitude of the regional models
    longitude_bounds: the maximum and minimum longitude of the regional models
    """
    # Define the 3 points on the west or east boundary. If west, we are on the 
    # minimum longitude, and if east we are on the maximum longitude.
    if boundary_name == "west":
        spherical_point_1 = np.array([np.max(radius_bounds), \
                                      np.max(latitude_bounds), \
                                      np.min(longitude_bounds)])
        spherical_point_2 = np.array([np.max(radius_bounds), \
                                      np.min(latitude_bounds), \
                                      np.min(longitude_bounds)])
        spherical_point_3 = np.array([np.min(radius_bounds), \
                                      np.max(latitude_bounds), \
                                      np.min(longitude_bounds)])

    elif boundary_name == "east":
        spherical_point_1 = np.array([np.max(radius_bounds), \
                                     np.max(latitude_bounds), \
                                     np.max(longitude_bounds)])
        spherical_point_2 = np.array([np.max(radius_bounds), \
                                     np.min(latitude_bounds), \
                                     np.max(longitude_bounds)])
        spherical_point_3 = np.array([np.min(radius_bounds), \
                                     np.max(latitude_bounds), \
                                     np.max(longitude_bounds)])

    else:
        raise Exception("Unknown boundary name: " + boundary_name)
    # Convert spherical points to Cartesian
    cartesian_point_1 = spherical_to_cartesian(spherical_point_1[0], \
                                               spherical_point_1[1], \
                                               spherical_point_1[2], \
                                               for_slice=False)
    cartesian_point_2 = spherical_to_cartesian(spherical_point_2[0], \
                                               spherical_point_2[1], \
                                               spherical_point_2[2], \
                                               for_slice=False)
    cartesian_point_3 = spherical_to_cartesian(spherical_point_3[0], \
                                               spherical_point_3[1], \
                                               spherical_point_3[2], \
                                               for_slice=False)
    # Calculate 2 in-plane orthogonal vectors using the 3 Cartesian points
    vector_1_2 = cartesian_point_2 - cartesian_point_1
    vector_1_3 = cartesian_point_3 - cartesian_point_1
    # Taking the cross product yields a vector normal to the model boundary
    normal_vector_to_plane = np.cross(vector_1_2, vector_1_3)
    # Normalize
    unit_normal = normal_vector_to_plane / np.linalg.norm(normal_vector_to_plane)

    return unit_normal    

####################################################################################################################################

""" Usage: This script requires 8 input arguments:
input_data: solution file for the global model (*.pvd)
output_directory: Where the .txt files for each boundary are saved
refinement_level: The number of mesh refinements in the global model
output_radius_resolution: The radial resolution of the regional slice (in meters)
output_lateral_resolution: The lateral resolution of the regional slice (in degrees)
radius_bounds: Array with the minimum and maximum radius (meters) of the regional model
latitude_bounds: Array with the minimum and maximum latitude (degrees) of the regional model
longitude_bounds: Array with the minimum and maximum longitude (degrees) of the regional model
"""

# Define the input arguments for the S2ORTS cookbook
input_data = "../../../cookbooks/initial-condition-S20RTS/output-S20RTS/solution.pvd"
output_directory = "./regional_velocity_files/"
refinement_level = 2
output_radius_resolution = 10e3 # 10 km radial resolution
output_lateral_resolution = 0.25 # 0.25 degree lateral resolution
radius_bounds = np.array([4770e3, 6360e3]) # Radius bounds of the regional model
latitude_bounds = np.array([-20, -55]) # Latitude bounds of the regional model
longitude_bounds = np.array([152, 210]) # Longitude bounds of the regional model

# Load in the global model
model = simple.OpenDataFile(input_data)

# Only load in the mesh-points and the velocity fields for computational efficiency, to load
# more fields, add entries to this array or just comment the line to load all solution fields.
model.PointArrays = ['Points', 'velocity']
model.UpdatePipeline() # Apply the filter

# Loop over the 4 lateral boundaries: west, east, north, south.
for boundary_name in np.array(["west", "east", "north", "south"]):
    # To prescribe the velocity on the boundaries of a 3D spherical chunk in ASPECT, the ASCII files must be 
    # named following this convention: chunk_3d_%s.%d.txt, where %s is the name of the model boundary, and 
    # %d is the timestep that the ASCII file will be applied. This python script is intended to be used for 
    # instantaneous models, and so %d is hardcoded to 0 when setting the variable output file name.
    output_datafile = output_directory + "chunk_3d_" + boundary_name + ".0.txt"
    
    # pvpython will save the data in a .csv file, which will not be formatted correctly for use in ASPECT.
    # This script creates a temporary file to store the pvpython output, then deletes the file later.
    temp_datafile = output_directory + boundary_name + "_boundary_paraview_slice_temp.csv"

    # For the east and west boundary, use the slice filter to extract velocities since these boundaries 
    # lie on great circles.
    if boundary_name == "west" or boundary_name == "east":
        cut_slice = simple.Slice(Input=model)
        cut_slice.SliceType.Normal = slice_plane_calculator(boundary_name, \
                                                            radius_bounds, \
                                                            latitude_bounds, \
                                                            longitude_bounds)
        cut_slice.UpdatePipeline()

    # For the north and south boundary, use the threshold filter at a constant latitude, since these
    # boundaries will not always lie on great circles.
    elif boundary_name == "north" or boundary_name == "south":
        # Create a calculator filter to determine the latitude in the global models
        pythonCalculator = simple.PythonCalculator(registrationName="pythonCalculator", Input=model)
        pythonCalculator.ArrayAssociation = "Point Data" # "Point Data" since the velocity is output on points

        # Copy all of the other variables into the calculator filter. If this is set to False, then the 
        # velocity field will not be passed through once the pythonCalculator is applied
        pythonCalculator.CopyArrays = True
        # Calculate the latitude
        pythonCalculator.Expression = "90 - np.arccos( points[:, 2] / (sqrt(points[:, 0]**2 + points[:, 1]**2 + points[:, 2]**2)) ) * 180 / np.pi"
        pythonCalculator.ArrayName = "latitude"
        pythonCalculator.UpdatePipeline() # Apply the filter

        # Create a threshold filter
        cut_slice = simple.Threshold(Input=pythonCalculator)
        cut_slice.Scalars = ("POINTS", "latitude") # Threshold the latitude variable
        cut_slice.ThresholdMethod = "Between" # Specify the 'between' method for the threshold filter

        # Threshold on either side of the maximum lat_bound (north) or the minimum lat_bound (south)
        # based on the refinement_level of the global models.
        if boundary_name == "north":
            cut_slice.LowerThreshold = np.max(latitude_bounds) - (180 / 2**(refinement_level + 1))
            cut_slice.UpperThreshold = np.max(latitude_bounds) + (180 / 2**(refinement_level + 1))
            
        if boundary_name == "south":
            cut_slice.LowerThreshold = np.min(latitude_bounds) - (180 / 2**(refinement_level + 1))
            cut_slice.UpperThreshold = np.min(latitude_bounds) + (180 / 2**(refinement_level + 1))
            
        cut_slice.UpdatePipeline()

    else:
        raise Exception("Unknown boundary name: " + boundary_name)

    # Create an array for the radius of the regional slice
    radius_array = np.arange(min(radius_bounds), max(radius_bounds) + output_radius_resolution, output_radius_resolution)

    # Create arrays for the latitude and longitude for the regional slices depending on which boundary is 
    # currently within the loop
    if boundary_name == "west":
        latitude_array = np.arange(min(latitude_bounds), max(latitude_bounds) + output_lateral_resolution, output_lateral_resolution)
        longitude_array = np.array([np.min(longitude_bounds)])
        
    elif boundary_name == "east":
        latitude_array = np.arange(min(latitude_bounds), max(latitude_bounds) + output_lateral_resolution, output_lateral_resolution)
        longitude_array = np.array([np.max(longitude_bounds)])
    
    elif boundary_name == "north":
        latitude_array = np.array([np.max(latitude_bounds)])
        longitude_array = np.arange(min(longitude_bounds), max(longitude_bounds) + output_lateral_resolution, output_lateral_resolution)
        
    elif boundary_name == "south":
        latitude_array = np.array([np.min(latitude_bounds)])
        longitude_array = np.arange(min(longitude_bounds), max(longitude_bounds) + output_lateral_resolution, output_lateral_resolution)

    # Calculate the cartesian coordinates on the regional slice
    cartesian_slice_coords = spherical_to_cartesian(radius_array, np.flip(latitude_array), longitude_array, for_slice=True)

    # Save Cartesian coordinates to .csv file for use later
    np.savetxt(temp_datafile, X=cartesian_slice_coords, delimiter=',', header="x,y,z", comments='')
    
    # This opens the .csv as a paraview object
    reader = simple.OpenDataFile(temp_datafile)
    reader.UpdatePipeline()

    # Create a TableToPoints filter that will be used later for storing interpolated data from the slices
    tabletopoints = servermanager.filters.TableToPoints()
    
    # Assign the TableToPoints object the points of the .csv file
    tabletopoints.XColumn = 'x'
    tabletopoints.YColumn = 'y'
    tabletopoints.ZColumn = 'z'
    tabletopoints.Input = reader
    tabletopoints.UpdatePipeline()

    # Take the cut_slice object and interpolate the fields onto the TableToPoints filter above.
    resample = simple.ResampleWithDataset(SourceDataArrays=cut_slice, DestinationMesh=tabletopoints)
    # Allow pvpython to compute the tolerance for resampling. If set to False, the user can provide
    # a tolerance with the line `resample.Tolerance`
    resample.ComputeTolerance = True
    resample.UpdatePipeline()

    # Take the points defined for the interpolated slice, and resave it in the .csv file.
    writer = simple.DataSetCSVWriter()
    writer.Input = resample
    writer.WriteTimeSteps = 1
    writer.FileName = temp_datafile      
    writer.UpdatePipeline()

    # Load the .csv file with the interpolated slice velocity data and extract both
    # the x, y, z coordinates and the x, y, and z components of the velocity.
    dataframe = pd.read_csv(temp_datafile)
    x_slice = dataframe["Points:0"].to_numpy()
    y_slice = dataframe["Points:1"].to_numpy()
    z_slice = dataframe["Points:2"].to_numpy()

    vx_slice = dataframe["velocity:0"].to_numpy()
    vy_slice = dataframe["velocity:1"].to_numpy()
    vz_slice = dataframe["velocity:2"].to_numpy()

    # Convert the cartesian points to spherical points
    r_slice, theta_slice, phi_slice = cartesian_to_spherical(x_slice, y_slice, z_slice)
    
    # Create the header and the arrays for saving into an ASCII file. Round to ensure that
    # each entry for the radius and latitude/longitude is consistent. The file will be named
    # and formatted in such a way that it can be directly applied to an ASPECT model.
    if boundary_name == "west" or boundary_name == "east":
        # ASPECT expects colatitude, so convert theta_slice to colatitude
        into_ASCII = np.array([np.round(r_slice, -3), \
                               np.deg2rad(np.round(90 - theta_slice, 2)), \
                               vx_slice, \
                               vy_slice, \
                               vz_slice]).T    
        header = "POINTS: " + str(len(radius_array)) + " " + str(len(latitude_array))

    elif boundary_name == "north" or boundary_name == "south":
        into_ASCII = np.array([np.round(r_slice, -3), \
                               np.deg2rad(np.round(phi_slice, 2)), \
                               vx_slice, \
                               vy_slice, \
                               vz_slice]).T 
        header = "POINTS: " + str(len(radius_array)) + " " + str(len(longitude_array))

    # Save into an ASCII file and remove the temporary .csv files
    os.remove(temp_datafile) 
    np.savetxt(fname=output_datafile, X=into_ASCII, fmt='%.5e', header=header)
    