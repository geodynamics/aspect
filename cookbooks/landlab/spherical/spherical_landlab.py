
from mpi4py import MPI
import json
import os
import numpy as np
import landlab
from landlab.components import LinearDiffuser
from landlab.components import FlowAccumulator
from landlab.components import StreamPowerEroder

from landlab.io.legacy_vtk import write_legacy_vtk

current_time = 0

comm = None

model_grid = None
elevation = None
linear_diffuser = None
flow_accumulator = None
stream_power_eroder = None

s2yr = 365.25 * 24 * 3600
timestep = 0
vtks = []

def cartesian_to_spherical(x, y, z):
    radius = np.sqrt(x**2 + y**2 + z**2)
    theta  = np.rad2deg(np.arccos(z / radius))
    phi    = np.rad2deg(np.arctan2(y, x))

    return radius, theta, phi

def initialize(comm_handle):
    if not comm_handle is None:
        # Convert the handle back to an MPI communicator
        global comm
        comm = MPI.Comm.f2py(comm_handle)

        rank = comm.Get_rank()
        size = comm.Get_size()

        print(f"Python: Hello from Rank {rank} of {size}")

        data = 1
        globalsum = comm.allreduce(data, op=MPI.SUM)
        if comm.rank == 0:
            print(f"\tPython: testing communication; sum {globalsum}")
    else:
        print("Python: running sequentially!")

def finalize():
    pass


# Run the Landlab simulation from the current time to end_time and return
# the new topographic elevation (in m) at each local node.
# dict_variable_name_to_value_in_nodes is a dictionary mapping variables
# (x velocity, y velocity, temperature, etc.) to an array of values in each
# node.
def update_until(end_time, dict_variable_name_to_value_in_nodes):
    global current_time, elevation, linear_diffuser, flow_accumulator, stream_power_eroder, timestep

    dt = end_time - current_time
    timestep += 1

    deposition_erosion = np.zeros(model_grid.number_of_nodes)

    x_velocity = dict_variable_name_to_value_in_nodes["x velocity"]
    y_velocity = dict_variable_name_to_value_in_nodes["y velocity"]
    z_velocity = dict_variable_name_to_value_in_nodes["z velocity"]

    radial_velocity = np.sqrt(x_velocity**2 + y_velocity**2 + z_velocity**2)

    # Substepping for surface processes
    if dt>0:
        n_substeps = 10
        sub_dt = dt / n_substeps
        for _ in range(n_substeps):

          # TODO:
          elevation_before = elevation.copy()

          #flow_accumulator.run_one_step()
          #stream_power_eroder.run_one_step(sub_dt)
          linear_diffuser.run_one_step(sub_dt)

          elevation += radial_velocity * sub_dt

          deposition_erosion += elevation - elevation_before
        pass

    filename = f"./output-spherical-landlab/landlab_{str(timestep).zfill(3)}.vtk"
    print("Writing output VTK file...", filename)
    vtk_file = write_legacy_vtk(path=filename, grid=model_grid, clobber=True, z_at_node=elevation)
    
    vtks.append((current_time, filename))

    if True:
        # write vtk.series file (ParaView supports legacy VTK in this format)
        with open("./output-spherical-landlab/landlab.vtk.series", "w") as f:
            series = {
                "file-series-version": "1.0",
                "files": [
                    {"name": os.path.basename(filename), "time": time / s2yr}
                    for time, filename in vtks
                ]
            }
            json.dump(series, f, indent=2)

    current_time = end_time
    print("Max elevation:", np.max(elevation), 
        "Min elevation:", np.min(elevation),
        "min deposition:", np.min(deposition_erosion),
        "max deposition:", np.max(deposition_erosion))
    
    return deposition_erosion.copy()

def set_mesh_information(dict_grid_information):
    global model_grid, elevation

    if not model_grid:
        print("* Creating Spherical Grid ...")
        #x_extent = dict_grid_information["Mesh X extent"]
        #y_extent = dict_grid_information["Mesh Y extent"]
        #spacing = dict_grid_information["Mesh Spacing"]

        model_grid = landlab.IcosphereGlobalGrid(radius=6371e3*0.99, mesh_densification_level=3)
        print(type(model_grid.x_of_node),type(model_grid.y_of_node),type(model_grid.z_of_node))
        for i in range(len(model_grid.x_of_node)):
            x=model_grid.x_of_node[i]
            y=model_grid.y_of_node[i]
            z=model_grid.z_of_node[i]
            print(i,x,y,z,np.sqrt(x*x+y*y+z*z))


        

        print("model grid nodes", len(model_grid.x_of_node))
        r, theta, phi = cartesian_to_spherical(model_grid.x_of_node, \
                                               model_grid.y_of_node, \
                                               model_grid.z_of_node)

        elevation = model_grid.add_zeros("topographic__elevation", at="node")
        elevation += 1e3 + np.random.rand(elevation.size) * 0.1e3
        #elevation[theta < 90] += 1000
        print(elevation)

        initialize_landlab_components(None)
        print("* Done")

# Return the x coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_x(grid_dictionary):
    global model_grid
    #dim = grid_dictionary["ASPECT Dimension"]
    #if dim == 2:
    #    print(np.unique(model_grid.x_of_node))
    #    return  np.unique(model_grid.x_of_node)
    #if dim == 3:
    return model_grid.x_of_node.copy()

# Return the y coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_y(grid_dictionary):
    global model_grid
    return model_grid.y_of_node.copy()
def get_grid_z(grid_dictionary):
    global model_grid
    return model_grid.z_of_node.copy()

def initialize_landlab_components(landlab_component_parameters):
    global model_grid, elevation, linear_diffuser, flow_accumulator, stream_power_eroder

    D    = 1e-4 #landlab_component_parameters["Hillslope coefficient"]
    K_sp = 1e-5 #landlab_component_parameters["Stream power erodibility coefficient"]
    flow_director = "Steepest" # Only flow director supported for HexGrids


    print("* Creating LinearDiffuser ... with D =", D)
    linear_diffuser = LinearDiffuser(model_grid, linear_diffusivity=D)
    print("* Creating FlowAccumulator ... with flow director =", flow_director)
    flow_accumulator = FlowAccumulator(model_grid, elevation, flow_director=flow_director)
    print("* Creating StreamPowerEroder ... with K_sp =", K_sp)
    stream_power_eroder = StreamPowerEroder(model_grid, K_sp=K_sp)

    pass
# Return the initial topography at the start of the simulation
# in each node.
def get_initial_topography(grid_dictionary):
    global elevation
    #dim = grid_dictionary["ASPECT Dimension"]
    #if dim == 2:
    #    return elevation[model_grid.y_of_node == 0]
    #if dim == 3:
    return elevation

