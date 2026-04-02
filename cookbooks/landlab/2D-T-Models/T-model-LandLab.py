
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

    slice_x_velocity = dict_variable_name_to_value_in_nodes["x velocity"]
    slice_y_velocity = dict_variable_name_to_value_in_nodes["y velocity"]

    # In the current setup, the y_velocity vector is ONLY nonzero along y=0. We need
    # to project the velocity outwards along all fixed values of y, assigning the same
    # velocity to all nodes that share the same x value. 

    # This code below first extract the velocities along y=0 (where they are correctly coming
    # from ASPECT), then it finds all of the unique y values in the LandLab mesh, and then assigned
    # the y=0 velocities to these y-values. This will only work for structured RasterGrids.

    vertical_velocity = np.zeros(model_grid.number_of_nodes)
    unique_x_values   = np.unique(model_grid.x_of_node)
    for x in unique_x_values:
        vertical_velocity[model_grid.x_of_node == x] = slice_y_velocity[unique_x_values == x]

    # Substepping for surface processes
    if dt>0:
        n_substeps = 10
        sub_dt = dt / n_substeps
        for _ in range(n_substeps):

          # TODO:
          elevation_before = elevation.copy()

          flow_accumulator.run_one_step()
          stream_power_eroder.run_one_step(sub_dt)
          linear_diffuser.run_one_step(sub_dt)

          elevation += vertical_velocity * sub_dt

          deposition_erosion += elevation - elevation_before
        pass

    filename = f"./output-2D-T-MODELS/landlab_{str(timestep).zfill(3)}.vtk"
    print("Writing output VTK file...", filename)
    vtk_file = write_legacy_vtk(path=filename, grid=model_grid, clobber=True, z_at_node=elevation)
    
    vtks.append((current_time, filename))

    if True:
        # write vtk.series file (ParaView supports legacy VTK in this format)
        with open("./output-2D-T-MODELS/landlab.vtk.series", "w") as f:
            series = {
                "file-series-version": "1.0",
                "files": [
                    {"name": os.path.basename(filename), "time": time / s2yr}
                    for time, filename in vtks
                ]
            }
            json.dump(series, f, indent=2)

    current_time = end_time
    print("Max elevation:", np.max(elevation), "Min elevation:", np.min(elevation))

    # This setup does not do any averaging of the topography across the LandLab mesh before returning
    # the change in elevation used to deform the ASPECT mesh. Averaging is straight forward to 
    # do going forward, and the commented line in the for loop below shows how this would be done.
    # However, for the purposes of this test I want to make sure the ASPECT mesh is
    # deforming to the exact same topography as the LandLab mesh where they intersect.  
    deposition_erosion_2d = np.zeros(len(np.unique(model_grid.x_of_node)))
    for x in unique_x_values:
        # deposition_erosion_2d[unique_x_values == x] = np.average(deposition_erosion[model_grid.x_of_node == x])
        deposition_erosion_2d = deposition_erosion[model_grid.y_of_node == 0]
    
    return deposition_erosion_2d

def set_mesh_information(dict_grid_information):
    global model_grid, elevation

    if not model_grid:
        print("* Creating RasterModelGrid ...")
        x_extent = 102.5e3
        y_extent = 100e3
        spacing  = 2500

        nrows = int(y_extent / spacing) + 1  # number of node rows
        ncols = int(x_extent / spacing) + 1  # number of node columns

        # For 2D models, shift the LandLab mesh so that it is centered on the y-axis.
        model_grid = landlab.RasterModelGrid((nrows, ncols), xy_spacing=(spacing, spacing), xy_of_lower_left=(0, -y_extent / 2))

        model_grid.set_closed_boundaries_at_grid_edges(right_is_closed=True, 
                                                       left_is_closed=False, 
                                                       top_is_closed=True, 
                                                       bottom_is_closed=True)

        print("* Creating topographic elevation ...")
        # Initialize topography array with zeros
        elevation = model_grid.add_zeros("topographic__elevation", at="node")
        # Assign 1000 m high topography with random noise
        np.random.seed(0)  # For reproducibility
        elevation += 1000 + np.random.rand(model_grid.number_of_nodes)

        initialize_landlab_components(None)
        print("* Done")

# Return the x coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_x(grid_dictionary):
    global model_grid
    return  np.unique(model_grid.x_of_node)

# Return the y coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_y(grid_dictionary):
    global model_grid
    return np.zeros(np.unique(model_grid.x_of_node).shape)

def initialize_landlab_components(grid_dictionary):
    global model_grid, elevation, linear_diffuser, flow_accumulator, stream_power_eroder

    D    = 1e-4 / s2yr
    K_sp = 1e-5 / s2yr
    flow_director = "D8" # Only flow director supported for HexGrids


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
    return elevation[model_grid.y_of_node == 0]
