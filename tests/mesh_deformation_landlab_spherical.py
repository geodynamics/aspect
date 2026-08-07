import numpy as np
import landlab
from landlab.components import LinearDiffuser

current_time    = 0
model_grid      = None
elevation       = None
linear_diffuser = None

def initialize(comm_handle):
    pass

def finalize():
    pass


# Run the Landlab simulation from the current time to end_time and return
# the new topographic elevation (in m) at each local node.
# dict_variable_name_to_value_in_nodes is a dictionary mapping variables
# (x velocity, y velocity, temperature, etc.) to an array of values in each
# node.
def update_until(aspect_solution_dict, aspect_auxiliary_dict):
    global current_time, linear_diffuser
    end_time = aspect_auxiliary_dict["ASPECT model time"]
    dt = end_time - current_time
    
    deposition_erosion = np.zeros(model_grid.number_of_nodes)

    print(f"Landlab running update_until: end_time = {end_time}", flush=True)

    x_velocity = aspect_solution_dict["x velocity"]
    y_velocity = aspect_solution_dict["y velocity"]
    z_velocity = aspect_solution_dict["z velocity"]

    if dt>0:
        n_substeps = 10
        sub_dt = dt / n_substeps
        for _ in range(n_substeps):
          elevation_before = elevation
          
          linear_diffuser.run_one_step(sub_dt)
          
          deposition_erosion += elevation - elevation_before
        pass
    
    current_time = end_time
    
    return deposition_erosion

def set_mesh_information(dict_grid_information):
    global model_grid, elevation, linear_diffuser

    if not model_grid:
        print("* Creating IcosphereGlobalGrid ...", flush=True)

        model_grid = landlab.IcosphereGlobalGrid(radius=20, mesh_densification_level=2)

        print("* The number of Landlab grid nodes is: ", len(model_grid.x_of_node), flush=True)
        print("* The node coordinates (in radians) are:", model_grid.phi_of_node, model_grid.theta_of_node, flush=True)

        print("* Creating topographic elevation ...", flush=True)
        elevation = model_grid.add_zeros("topographic__elevation", at="node")
        elevation[model_grid.theta_of_node < np.pi/2] += 1.0

        D = 0.01
        print("* Creating LinearDiffuser ... with D =", D, flush=True)
        linear_diffuser = LinearDiffuser(model_grid, linear_diffusivity=D)

        print("* Done initializing Landlab mesh.", flush=True)

# Return the x coordinates of the locally owned nodes on this
# MPI rank.
def get_grid_x(grid_dictionary):
    global model_grid
    return model_grid.x_of_node.copy()

# Return the y coordinates of the locally owned nodes on this
# MPI rank.
def get_grid_y(grid_dictionary):
    global model_grid
    return model_grid.y_of_node.copy()

# Return the z coordinates of the locally owned nodes on this
# MPI rank.
def get_grid_z(grid_dictionary):
    global model_grid
    return model_grid.z_of_node.copy()

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
    print("grid coordinates:", get_grid_x(0), get_grid_y(0), flush=True)

    dt = 0.1
    for n in range(3):
        data = {}
        data["x velocity"] = np.zeros(model_grid.number_of_nodes)
        data["y velocity"] = np.zeros(model_grid.number_of_nodes)
        data["z velocity"] = np.zeros(model_grid.number_of_nodes)

        time = {}
        time["ASPECT model time"] = n*dt
        update_until(data, time)
        write_output()
