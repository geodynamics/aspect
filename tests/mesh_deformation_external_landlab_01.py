print("mesh_deformation_external_landlab_01.py")

from mpi4py import MPI
import numpy as np
import landlab
from landlab.components import LinearDiffuser

current_time = 0

comm = None

mg = None
elevation = None
linear_diffuser = None

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

def finalize():
    pass


# Run the Landlab simulation from the current time to end_time and return
# the new topographic elevation (in m) at each local node.
# dict_variable_name_to_value_in_nodes is a dictionary mapping variables
# (x velocity, y velocity, temperature, etc.) to an array of values in each
# node.
def update_until(end_time, dict_variable_name_to_value_in_nodes):
    global current_time
    dt = end_time - current_time
    
    deposition_erosion = np.zeros(mg.number_of_nodes)

    print(f"update_until: end_time = {end_time}")

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
          elevation_before = elevation
          
          linear_diffuser.run_one_step(sub_dt)
          
          deposition_erosion += elevation - elevation_before
        pass
    
    current_time = end_time
    
    return deposition_erosion

def set_mesh_information(dict_grid_information):
    global mg, elevation, linear_diffuser

    if not mg:
        print("* Creating HexModelGrid ...")
        mg = landlab.HexModelGrid((5, 6), node_layout="rect", spacing=0.2)
        print("* Creating topographic elevation ...")
        elevation = mg.add_zeros("topographic__elevation", at="node")
        # Add reproducible noise to topography:
        np.random.seed(42)
        elevation += 1000 + np.random.rand(elevation.size) / 10.0
        print("\tnumber of nodes:", mg.number_of_nodes)
        print("\tnode coordinates:", mg.node_x, mg.node_y)

        D = 0.01 # m2
        print("* Creating LinearDiffuser ... with D =", D)
        linear_diffuser = LinearDiffuser(mg, linear_diffusivity=D)

        print("* Done")

# Return the x coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_x(grid_id):
    global mg
    return mg.node_x

# Return the y coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_y(grid_id):
    global mg
    return mg.node_y

# Return the initial topography at the start of the simulation
# in each node.
def get_initial_topography(grid_id):
    global elevation
    return elevation


def write_output():
    pass    

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    initialize(MPI.Comm.py2f(comm))

    set_mesh_information({})
    print("grid coordinates:", get_grid_x(0), get_grid_y(0))

    dt = 0.1
    for n in range(3):
        data = {}
        data["x velocity"] = np.zeros(mg.number_of_nodes)
        data["y velocity"] = np.zeros(mg.number_of_nodes)
        data["z velocity"] = np.zeros(mg.number_of_nodes)
        update_until(n*dt, data)
        write_output()
