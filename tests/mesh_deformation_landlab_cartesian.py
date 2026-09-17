import numpy as np
import landlab
from landlab.components import LinearDiffuser

from landlab_template import LandlabTemplate


class MyAspectLandlabModel(LandlabTemplate):

    # Run the Landlab simulation from the current time to end_time and return
    # the new topographic elevation (in m) at each local node.
    # dict_variable_name_to_value_in_nodes is a dictionary mapping variables
    # (x velocity, y velocity, temperature, etc.) to an array of values in each
    # node.
    def update_until(self, ASPECT_solution_at_Landlab_nodes_dict, ASPECT_additional_info_dict):
        end_time = ASPECT_additional_info_dict["ASPECT model time"]
        dt = end_time - self.current_time
        
        deposition_erosion = np.zeros(self.model_grid.number_of_nodes)

        print(f"Landlab running update_until: end_time = {end_time}", flush=True)

        x_velocity = ASPECT_solution_at_Landlab_nodes_dict["x velocity"]
        y_velocity = ASPECT_solution_at_Landlab_nodes_dict["y velocity"]
        z_velocity = ASPECT_solution_at_Landlab_nodes_dict["z velocity"]

        if dt>0:
            n_substeps = 10
            sub_dt = dt / n_substeps
            for _ in range(n_substeps):
                elevation_before = self.elevation.copy()
                
                self.linear_diffuser.run_one_step(sub_dt)
                
                deposition_erosion += self.elevation - elevation_before
        
        self.current_time = end_time
        
        return deposition_erosion

    def set_mesh_information(self, grid_dictionary):

        if not self.model_grid:
            print("* Creating HexModelGrid ...", flush=True)
            self.model_grid = landlab.HexModelGrid((5, 6), node_layout="rect", spacing=0.2)

            print("* The number of Landlab grid nodes is: ", len(self.model_grid.x_of_node), flush=True)
            print("* The node coordinates (in meters) are:", self.model_grid.node_x, self.model_grid.node_y, flush=True)

            print("* Creating topographic elevation ...", flush=True)
            self.elevation = self.model_grid.add_zeros("topographic__elevation", at="node")

            # Add reproducible noise to topography:
            np.random.seed(42)
            self.elevation += 1000 + np.random.rand(self.elevation.size) / 10.0

            self.initialize_landlab_components()


    def initialize_landlab_components(self):
        D = 0.01 # m2
        print("* Creating LinearDiffuser ... with D =", D, flush=True)
        self.linear_diffuser = LinearDiffuser(self.model_grid, linear_diffusivity=D)

        print("* Done initializing Landlab mesh.", flush=True)


model = MyAspectLandlabModel()
model.export_ASPECT_functions(model, globals())

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    model.initialize(MPI.Comm.py2f(comm))

    model.set_mesh_information({})
    print("grid coordinates:", model.get_grid_x(0), model.get_grid_y(0))

    dt = 0.1
    for n in range(3):
        data = {}
        data["x velocity"] = np.zeros(model.model_grid.number_of_nodes)
        data["y velocity"] = np.zeros(model.model_grid.number_of_nodes)
        data["z velocity"] = np.zeros(model.model_grid.number_of_nodes)

        time = {}
        time["ASPECT model time"] = n*dt
        model.update_until(data, time)
        model.write_output()
