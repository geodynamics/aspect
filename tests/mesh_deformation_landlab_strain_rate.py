import numpy as np
import landlab

current_time    = 0
model_grid      = None

def initialize(comm_handle):
    pass

def finalize():
    pass


# Check to see that the strain rate is being sent to Landlab from ASPECT.
def update_until(aspect_solution_dict, aspect_auxiliary_dict):
    global current_time
    end_time = aspect_auxiliary_dict["ASPECT model time"]
    dt = end_time - current_time
    current_time = end_time

    deposition_erosion = np.zeros(model_grid.number_of_nodes)

    strain_rate = aspect_solution_dict["strain rate"]
    if dt > 0:
        print(f"The Strain Rate in Landlab is: {strain_rate}", flush=True)
    
    return deposition_erosion

def set_mesh_information(dict_grid_information):
    global model_grid

    if not model_grid:
        x_extent = 80e3
        y_extent = 80e3
        spacing  = 10e3

        nrows = int(y_extent / spacing) + 1 # number of node rows
        ncols = int(x_extent / spacing) + 1 # number of node columns

        model_grid = landlab.RasterModelGrid((nrows, ncols), xy_spacing=(spacing, spacing), xy_of_lower_left=(0, 0))

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
    global model_grid
    return np.zeros(model_grid.number_of_nodes)


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
        data["x velocity"]  = np.zeros(model_grid.number_of_nodes)
        data["y velocity"]  = np.zeros(model_grid.number_of_nodes)
        data["z velocity"]  = np.zeros(model_grid.number_of_nodes)
        data["strain rate"] = np.zeros(model_grid.number_of_nodes)

        time = {}
        time["ASPECT model time"] = n*dt
        update_until(data, time)
        write_output()
