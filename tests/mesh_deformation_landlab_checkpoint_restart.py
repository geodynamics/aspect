import numpy as np
import landlab
import os
from landlab.components import LinearDiffuser
from landlab.io.native_landlab import save_grid, load_grid
from landlab.io.legacy_vtk import write_legacy_vtk

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
    global current_time, linear_diffuser, model_grid, elevation
    end_time = aspect_auxiliary_dict["ASPECT model time"]
    dt = end_time - current_time
    
    deposition_erosion = np.zeros(model_grid.number_of_nodes)

    print(f"Landlab running update_until: end_time = {end_time}", flush=True)

    x_velocity = aspect_solution_dict["x velocity"]
    y_velocity = aspect_solution_dict["y velocity"]

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
        print("* Creating HexModelGrid ...", flush=True)
        x_extent = 100e3
        y_extent = 100e3
        spacing  = 10e3
        nrows = x_extent / spacing + 1
        ncols = y_extent / spacing + 1
        model_grid = landlab.RasterModelGrid((int(nrows), int(ncols)), xy_spacing=(spacing, spacing))

        print("* The number of Landlab grid nodes is: ", len(model_grid.x_of_node), flush=True)
        print("* The node coordinates (in meters) are:", model_grid.node_x, model_grid.node_y, flush=True)

        print("* Creating topographic elevation ...", flush=True)
        elevation = model_grid.add_zeros("topographic__elevation", at="node")

        # Add reproducible noise to topography:
        np.random.seed(42)
        elevation += 1000 + np.random.rand(elevation.size) / 10.0

        print("* Done initializing Landlab mesh.", flush=True)

        # Initialize Landlab components after creating the mesh and topography
        initialize_landlab_components()


def initialize_landlab_components():
    global linear_diffuser, model_grid, elevation
    D = 1e-2 / 60 / 60 / 24 / 365.25 # m2
    print("* Creating LinearDiffuser ... with D =", D, flush=True)
    linear_diffuser = LinearDiffuser(model_grid, linear_diffusivity=D)

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


def checkpoint_model_grid(checkpoint_dict):
    """
    Checkpoint the Landlab model grid by saving it to a file. This function is called when
    checkpointing the ASPECT model.
    """
    global model_grid

    # Extract checkpointing information from the checkpoint dictionary
    checkpoint_index      = checkpoint_dict["Current checkpoint ID"]
    output_directory      = checkpoint_dict["Output directory"]

    # Create LandLab checkpoint directory within the ASPECT output directory.
    os.makedirs(output_directory, exist_ok=True)

    print(f"Checkpointing the LandLab model grid with index {checkpoint_index} ...", flush=True)
    filename = os.path.join(output_directory, f"landlab_checkpoint_{str(checkpoint_index).zfill(2)}.grid")
    save_grid(model_grid, filename, clobber=True)
    pass

def load_model_grid(checkpoint_dict):
    """
    Load the Landlab model grid from a file. This function is called when
    restarting the ASPECT model from a checkpoint.
    """
    global model_grid, elevation

    restart_checkpoint_id = checkpoint_dict["Resume checkpoint ID"]
    print(f"Restarting from checkpoint {restart_checkpoint_id} ...", flush=True)

    # Extract checkpointing information from the checkpoint dictionary
    output_directory = checkpoint_dict["Output directory"]

    # Load the LandLab grid from the checkpoint file corresponding to the checkpoint index.
    print(f"Loading the LandLab model grid from checkpoint {restart_checkpoint_id} ...", flush=True)
    filename = os.path.join(output_directory, f"landlab_checkpoint_{str(restart_checkpoint_id).zfill(2)}.grid")

    model_grid = load_grid(filename)
    elevation = model_grid.at_node["topographic__elevation"]

    # We need to initialize the components after loading the LandLab grid, since these
    # are not stored.
    initialize_landlab_components()
    pass

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
