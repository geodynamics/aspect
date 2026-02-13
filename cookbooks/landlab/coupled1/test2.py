print("Hello from test2.py")

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

uplift_rate = 1e-3 # m/yr
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

    x_velocity = dict_variable_name_to_value_in_nodes["x velocity"]
    y_velocity = dict_variable_name_to_value_in_nodes["y velocity"]
    z_velocity = dict_variable_name_to_value_in_nodes["z velocity"]

    if dt>0:
        n_substeps = 10
        sub_dt = dt / n_substeps
        for _ in range(n_substeps):
          
          # TODO:
          #uplift(z_velocity * sub_dt)
          #advect(x_velocity * sub_dt, y_velocity * sub_dt)
          elevation_before = elevation.copy()
          
          # Run the landlab components: First route the water, use stream power to erode based on
          # the routing of the water, then diffuse the topography
          flow_accumulator.run_one_step()
          stream_power_eroder.run_one_step(sub_dt)
          linear_diffuser.run_one_step(sub_dt)

          # Uplift the topography at a constant rate.
          elevation[model_grid.core_nodes] += uplift_rate * sub_dt
          
          deposition_erosion += elevation - elevation_before
        pass
    
    current_time = end_time

    filename = f"./output/landlab_{str(timestep).zfill(3)}.vtk"
    print("Writing output VTK file...", filename)
    vtk_file = write_legacy_vtk(path=filename, grid=model_grid, clobber=True)
    vtks.append((current_time, filename))

    if True:
        # write vtk.series file (ParaView supports legacy VTK in this format)
        with open("./output/landlab.vtk.series", "w") as f:
            series = {
                "file-series-version": "1.0",
                "files": [
                    {"name": os.path.basename(filename), "time": time}
                    for time, filename in vtks
                ]
            }
            json.dump(series, f, indent=2)
    
    print ("deposition/erosion:", np.linalg.norm(deposition_erosion))
    return deposition_erosion
    # return elevation

def set_mesh_information(dict_grid_information):
    global model_grid, elevation, linear_diffuser, flow_accumulator, stream_power_eroder

    if not model_grid:
        print("* Creating HexModelGrid ...")
        x_extent = 100e3 # 100 km 
        y_extent = 100e3 # 100 km 

        spacing = 1000 # 1 km spacing
        nrows = int(y_extent / spacing)  # number of node rows
        ncols = int(x_extent / spacing)  # number of node columns

        model_grid = landlab.HexModelGrid((nrows, ncols), spacing=spacing)

        print("* Creating topographic elevation ...")
        # Initialize topography array with zeros
        elevation = model_grid.add_zeros("topographic__elevation", at="node")
        # Assign 1000 m high topography with random noise
        elevation += 1000 + np.random.rand(elevation.size) / 1000.0
        print("\tnumber of nodes:", model_grid.number_of_nodes)

        D = 0.01 # m2
        K_sp = 1e-4 # Stream power erodibility coefficient
        flow_director = "Steepest" # Only flow director supported for HexGrids
        print("* Creating LinearDiffuser ... with D =", D)
        linear_diffuser = LinearDiffuser(model_grid, linear_diffusivity=D)
        print("* Creating FlowAccumulator ... with flow director =", flow_director)
        flow_accumulator = FlowAccumulator(model_grid, elevation, flow_director=flow_director)
        print("* Creating StreamPowerEroder ... with K_sp =", K_sp)
        stream_power_eroder = StreamPowerEroder(model_grid, K_sp=K_sp)

        print("* Done")

# Return the x coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_x(grid_id):
    global model_grid
    return model_grid.node_x

# Return the y coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_y(grid_id):
    global model_grid
    return model_grid.node_y

# Return the initial topography at the start of the simulation
# in each node.
def get_initial_topography(grid_id):
    global elevation
    return elevation


def write_output(timestep):
    # Write the grid to vtk
    print("Writing output VTK file...")
    vtk_file = write_legacy_vtk(path=f"./output_{str(timestep).zfill(3)}.vtk", grid=model_grid, clobber=True)

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    initialize(MPI.Comm.py2f(comm))

    set_mesh_information({})

    dt = 10.0e3 # 100 kyr timesteps, run for 1 Myr
    for n in range(10):
        data = {}
        data["x velocity"] = np.zeros(model_grid.number_of_nodes)
        data["y velocity"] = np.zeros(model_grid.number_of_nodes)
        data["z velocity"] = np.zeros(model_grid.number_of_nodes)
        update_until(n*dt, data)
        write_output(n)
