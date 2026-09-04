import numpy as np
import landlab
from landlab.components import LinearDiffuser
from landlab.components import FlowAccumulator, ErosionDeposition

current_time    = 0
model_grid      = None
elevation       = None
linear_diffuser = None
flow_accumulator = None
erosion_deposition = None
s2yr = 60 * 60 * 24 * 365 * 2  # seconds in 2 years

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
    global current_time, linear_diffuser, flow_accumulator, erosion_deposition, model_grid, elevation
    end_time = aspect_auxiliary_dict["ASPECT model time"]
    dt = end_time - current_time
    
    deposition_erosion = np.zeros(model_grid.number_of_nodes)

    print(f"Landlab running update_until: end_time = {end_time}", flush=True)

    slice_y_velocity = aspect_solution_dict["y velocity"]

    vertical_velocity = np.zeros(model_grid.number_of_nodes)
    unique_x_values = np.unique(model_grid.x_of_node)
    for x in unique_x_values:
        vertical_velocity[model_grid.x_of_node == x] = slice_y_velocity[unique_x_values == x]

    if dt>0:
        n_substeps = 10
        sub_dt = dt / n_substeps
        for _ in range(n_substeps):
          elevation_before = elevation.copy()

          elevation += vertical_velocity * sub_dt

          flow_accumulator.run_one_step()
          erosion_deposition.run_one_step(sub_dt)
          
          linear_diffuser.run_one_step(sub_dt)
          
          deposition_erosion += elevation - elevation_before
        pass
    
    current_time = end_time
    
    deposition_erosion_2d = np.zeros(len(np.unique(model_grid.x_of_node)))
    unique_x_values = np.unique(model_grid.x_of_node)

    for x in unique_x_values:
        deposition_erosion_2d[unique_x_values == x] = np.average(deposition_erosion[model_grid.x_of_node == x])
    return deposition_erosion_2d

def set_mesh_information(dict_grid_information):
    global model_grid, elevation, linear_diffuser, flow_accumulator, erosion_deposition, elevation

    if not model_grid:
        print("* Creating RasterModelGrid ...", flush=True)

        x_extent = 200e3
        y_extent = 100e3

        spacing = 2.5e3

        nrows = int(y_extent / spacing) + 1
        ncols = int(x_extent / spacing) + 1
        model_grid = landlab.RasterModelGrid((nrows, ncols), xy_spacing=(spacing,spacing), xy_of_lower_left=(0.0, -y_extent / 2))

        print("* The number of Landlab grid nodes is: ", len(model_grid.x_of_node), flush=True)
        print("* The node coordinates (in meters) are:", model_grid.node_x, model_grid.node_y, flush=True)

        print("* Creating topographic elevation ...", flush=True)
        elevation = model_grid.add_zeros("topographic__elevation", at="node")

        # Add reproducible noise to topography:
        np.random.seed(42)
        elevation += np.random.rand(elevation.size) * 5.0

        D = 0.01 # m2
        print("* Creating LinearDiffuser ... with D =", D, flush=True)
        linear_diffuser = LinearDiffuser(model_grid, linear_diffusivity=D / s2yr)
        flow_accumulator = FlowAccumulator(model_grid)
        erosion_deposition = ErosionDeposition(model_grid, m_sp=0.4, n_sp=1.0, K=1e-6 / s2yr)

        print("* Done initializing Landlab mesh.", flush=True)

# Return the x coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_x(grid_id):
    global model_grid
    return np.unique(model_grid.x_of_node)

# Return the y coordinates of the locally owned nodes on this
# MPI rank. grid_id is always 0.
def get_grid_y(grid_id):
    global model_grid
    return np.zeros(np.unique(model_grid.x_of_node).shape)

# Return the initial topography at the start of the simulation
# in each node.
def get_initial_topography(grid_id):
    global elevation
    return elevation[model_grid.y_of_node == 0]


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
        data["x velocity"] = np.zeros(model_grid.number_of_nodes)
        data["y velocity"] = np.zeros(model_grid.number_of_nodes)
        data["z velocity"] = np.zeros(model_grid.number_of_nodes)

        time = {}
        time["ASPECT model time"] = n*dt
        update_until(data, time)
        write_output()
