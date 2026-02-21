print("mesh_deformation_external_landlab_02.py")

from mpi4py import MPI
import numpy as np

current_time = 0

comm = None


px = None
py = None
elevation = None

def initialize(comm_handle):
    global px, py, elevation

    if not comm_handle is None:
        # Convert the handle back to an MPI communicator
        global comm
        comm = MPI.Comm.f2py(comm_handle)

        rank = comm.Get_rank()
        size = comm.Get_size()


        data = 1
        globalsum = comm.allreduce(data, op=MPI.SUM)
        if comm.rank == 0:
            print(f"Python: comm size = {size}, communication test = {globalsum}")
    
        if comm.rank == 0:
            px = np.array([0.1, 0.3, 0.4])
            py = np.array([0.2, 0.4, 0.5])

        if comm.rank == 1:
            px = np.array([0.5, 0.7])
            py = np.array([0.6, 0.8])

        elevation = np.zeros(px.size)

def finalize():
    pass


# Run the Landlab simulation from the current time to end_time and return
# the new topographic elevation (in m) at each local node.
# dict_variable_name_to_value_in_nodes is a dictionary mapping variables
# (x velocity, y velocity, temperature, etc.) to an array of values in each
# node.
def update_until(end_time, dict_variable_name_to_value_in_nodes):

    print(f"update_until: end_time = {end_time}")
   
    return elevation

def set_mesh_information(dict_grid_information):
    pass

# Return the x coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_x(grid_id):
    global px
    return px

# Return the y coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_y(grid_id):
    global py
    return py

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
        data["x velocity"] = np.zeros(px.size)
        data["y velocity"] = np.zeros(px.size)
        data["z velocity"] = np.zeros(px.size)
        update_until(n*dt, data)
        write_output()
