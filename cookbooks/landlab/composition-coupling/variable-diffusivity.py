
from mpi4py import MPI
import json
import os
import numpy as np
import landlab
from landlab.components import LinearDiffuser

from landlab.io.legacy_vtk import write_legacy_vtk

current_time = 0

comm = None

model_grid = None
elevation = None
linear_diffuser = None
Diffusivity = None

s2yr = 60 * 60 * 24 * 365.25
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
    global current_time, elevation, linear_diffuser, timestep, Diffusivity

    dt = end_time - current_time
    timestep += 1

    deposition_erosion = np.zeros(model_grid.number_of_nodes)

    # Extract the velocity along the ASPECT model, which is located at y=0 on the landlab mesh.
    slice_x_velocity = dict_variable_name_to_value_in_nodes["x velocity"]
    slice_y_velocity = dict_variable_name_to_value_in_nodes["y velocity"]
    # slice_z_velocity = dict_variable_name_to_value_in_nodes["z velocity"]

    # Likewise, extract the composition along the ASPECT surface. These exactly match the names
    # of the compositions provided to ASPECT.
    slice_weak_composition   = dict_variable_name_to_value_in_nodes["weak"]
    slice_strong_composition = dict_variable_name_to_value_in_nodes["strong"]

    # Create empty arrays to project the velocity and composition values out from y=0 to all
    # nodes on the landlab mesh.
    vertical_velocity = np.zeros(model_grid.number_of_nodes)
    strong_composition = np.zeros(model_grid.number_of_nodes)
    weak_composition   = np.zeros(model_grid.number_of_nodes)

    unique_x_values   = np.unique(model_grid.x_of_node)
    for x in unique_x_values:
        vertical_velocity[model_grid.x_of_node == x]  = slice_y_velocity[unique_x_values == x]
        strong_composition[model_grid.x_of_node == x] = slice_strong_composition[unique_x_values == x]
        weak_composition[model_grid.x_of_node == x]   = slice_weak_composition[unique_x_values == x]

    # Modify the diffusivity based on the composition, making it very high for "weak" and low for "strong".
    Diffusivity[strong_composition >= 0.5] = 1e-10 / s2yr
    Diffusivity[weak_composition >= 0.5]   = 10 / s2yr

    # Substepping for surface processes
    if dt>0:
        n_substeps = 10
        sub_dt = dt / n_substeps
        for _ in range(n_substeps):

          elevation_before = elevation.copy()

          linear_diffuser.run_one_step(sub_dt)

          elevation += vertical_velocity * sub_dt

          deposition_erosion += elevation - elevation_before
        pass

    current_time = end_time
    print("Max elevation:", np.max(elevation), "Min elevation:", np.min(elevation))

    # Return the change in the topography along y=0, where the ASPECT mesh is located.
    deposition_erosion_2d = np.zeros(len(np.unique(model_grid.x_of_node)))
    for x in unique_x_values:
        deposition_erosion_2d = deposition_erosion[model_grid.y_of_node == 0]
    
    return deposition_erosion_2d

def set_mesh_information(dict_grid_information):
    global model_grid, elevation

    if not model_grid:
        print("* Creating RasterModelGrid ...")
        x_extent = 100e3
        y_extent = 100e3
        spacing  = 500.0

        nrows = int(y_extent / spacing) + 1 # number of node rows
        ncols = int(x_extent / spacing) + 1 # number of node columns

        model_grid = landlab.RasterModelGrid((nrows, ncols), xy_spacing=(spacing, spacing), xy_of_lower_left=(0, -y_extent / 2))

        print("* Creating topographic elevation ...")
        # Initialize topography array with zeros
        elevation = model_grid.add_zeros("topographic__elevation", at="node")

        # Add a large 20 km high triangular mountain in the middle of the LandLab domain.
        topo_height = 20e3
        left_x_arr = np.array([25e3, 50e3])
        left_y_arr = np.array([0.0, topo_height])

        right_x_arr = np.array([50e3, 75e3])
        right_y_arr = np.array([topo_height, 0.0])

        left_m, left_b = np.polyfit(left_x_arr, left_y_arr, deg=1)
        right_m, right_b = np.polyfit(right_x_arr, right_y_arr, deg=1)

        elevation[model_grid.x_of_node <= 50e3] = left_m * model_grid.x_of_node[model_grid.x_of_node <= 50e3] + left_b
        elevation[model_grid.x_of_node > 50e3]  = right_m * model_grid.x_of_node[model_grid.x_of_node > 50e3] + right_b

        elevation[model_grid.x_of_node < np.min(left_x_arr)] = 0.0
        elevation[model_grid.x_of_node > np.max(right_x_arr)] = 0.0

        # Close all boundaries
        model_grid.set_closed_boundaries_at_grid_edges(right_is_closed=True, 
                                                       left_is_closed=True, 
                                                       top_is_closed=True, 
                                                       bottom_is_closed=True)
                                                       
        print("\tnumber of nodes:", model_grid.number_of_nodes)
        initialize_landlab_components(None)

        print("* Done")

# If the ASPECT mesh is 2D, then we only want to return the unique values of x on the LandLab mesh,
# the logic here is probably different if the LandLab mesh is not a raster.
def get_grid_x(grid_dictionary):
    global model_grid
    return np.unique(model_grid.x_of_node)

# If the ASPECT mesh is 2D, then we only return an array of 0s for the y values, since this is where
# the ASPECT surface is located.
def get_grid_y(grid_dictionary):
    global model_grid
    return np.zeros(np.unique(model_grid.x_of_node).shape)

def initialize_landlab_components(landlab_component_parameters):
    global model_grid, elevation, linear_diffuser, flow_accumulator, stream_power_eroder, Diffusivity

    D = 1e-10 / s2yr
    
    # If the diffusivity is mapped to an array with the same length as the number of nodes on the landlab
    # mesh, landlab will use those values as the diffusivity at each node. If these values are updated,
    # landlab will also track this.
    Diffusivity = model_grid.add_zeros("linear_diffusivity", at="node")
    Diffusivity += D

    # print("* Creating LinearDiffuser ... with D =", D)
    linear_diffuser = LinearDiffuser(model_grid, linear_diffusivity=Diffusivity)
    pass

# Return the initial topography along y=0, where the ASPECT surface is located.
def get_initial_topography(grid_dictionary):
    global elevation
    return elevation[model_grid.y_of_node == 0]


def write_output(timestep):
    # Write the grid to vtk
    print("Writing output VTK file...")
    vtk_file = write_legacy_vtk(path=f"./output_{str(timestep).zfill(3)}.vtk", grid=model_grid, clobber=True)

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    initialize(MPI.Comm.py2f(comm))

    set_mesh_information({})
    initialize_landlab_components(None)

     # Main simulation loop
    for n in range(100):
        data = {}
        data["x velocity"] = np.zeros(model_grid.number_of_nodes)
        data["y velocity"] = np.zeros(model_grid.number_of_nodes)
        data["z velocity"] = np.zeros(model_grid.number_of_nodes)
        update_until(n*dt, data)
        write_output(n)
