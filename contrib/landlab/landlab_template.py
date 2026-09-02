"""
Template class used for creating a LandLab script that is then used in the ASPECT-Landlab 
coupling. This can be imported and the functions can be modified to tailor the needs of the 
user for productions models or simple tests. 
"""

import inspect
import json
import os
import numpy as np
import landlab
from landlab.io.native_landlab import save_grid, load_grid
from landlab.io.legacy_vtk import write_legacy_vtk


class LandlabTemplate:
    """Template class for LandLab scripts used in ASPECT-Landlab coupling."""

    # Below is a list of functions that are explicitly called by ASPECT. Because
    # ASPECT will always call them with an expected argument signature, we check
    # that these functions retain their signature when overridden in a subclass.
    function_explicitly_called_by_aspect = (
        "initialize",
        "set_mesh_information",
        "update_until",
        "checkpoint_model_grid",
        "load_model_grid",
        "get_initial_topography",
        "write_output",
        "get_grid_x",
        "get_grid_y",
        "get_grid_z",
    )

    def __init_subclass__(sub_class, **kwargs):
        super().__init_subclass__(**kwargs)

        for function_name in LandlabTemplate.function_explicitly_called_by_aspect:
            # Only check functions that the subclass is directly overriding. If functions
            # are not overridden, the base class implementation is used and the
            # signature will be correct.
            if function_name not in sub_class.__dict__:
                continue

            # Get the base and subclass functions to compare their signatures.
            base_class_method = getattr(LandlabTemplate, function_name)
            sub_class_method = sub_class.__dict__[function_name]

            # Extract the signatures of the base and subclass functions.
            base_class_signature = inspect.signature(base_class_method)
            sub_class_signature = inspect.signature(sub_class_method)

            # Check if they are identical.
            if sub_class_signature != base_class_signature:
                raise TypeError(
                    f"Invalid override for '{function_name}' in class "
                    f"'{sub_class.__name__}': expected signature {base_class_signature}, "
                    f"got {sub_class_signature}."
                )

    def __init__(self):
        self.current_time = 0.0
        self.comm = None

        self.model_grid          = None
        self.elevation           = None
        self.horizontal_velocity = None

        # Use the same second to year conversion as in ASPECT
        self.seconds_in_year = 60 * 60 * 24 * 365.2425
        self.timestep = 0
        self.vtks = []

    def initialize(self, comm_handle):
        """Called once by ASPECT at startup."""
        if comm_handle is not None:
            self.comm = MPI.Comm.f2py(comm_handle)
            rank = self.comm.Get_rank()
            size = self.comm.Get_size()
            print(f"Python: Hello from Rank {rank} of {size}")

            globalsum = self.comm.allreduce(1, op=MPI.SUM)
            if self.comm.rank == 0:
                print(f"\tPython: testing communication; sum {globalsum}")
        else:
            print("Python: running sequentially!")

    # ---------------------------------------------------------------------------
    # The following functions will need to be modified based on the needs of the user:
    # - set_mesh_information: create the model grid and initialize the elevation field.
    #   Here the user can define the geometry and resolution of the landlab mesh.
    #
    # - initialize_landlab_components: initialize any LandLab components. Here the user
    #   chooses which Landlab components will influence the evolution of the topograhy, and
    #   initializes them with the necessary parameters.
    #
    # - update_until: advance the LandLab model through time based on the ASPECT timestep.
    #   Here the user needs to define how the Landlab components influence the evolution of
    #   the topography, as well as how fields from ASPECT (composition, temperature etc.,)
    #   influence the Landlab model.
    # ---------------------------------------------------------------------------
    def set_mesh_information(self, grid_dictionary):
        """
        Create the Landlab model grid, initialize the elevation field, and any other grid fields needed for
        the Landlab model. Call the function initialize_landlab_components() at the end of this function to 
        initialize the Landlab components after the grid and fields have been created.
        ** Derived classes must override this function in custom scripts. **
        """
        raise NotImplementedError("LandlabTemplate.set_mesh_information() must be overridden in your custom script. " \
                                  "See the doc string in the template file for information on the recommended " \
                                  "way to write this function.")



    def initialize_landlab_components(self):
        """
        Initialize the Landlab components that will be used to evolve the topography over time.
        ** Derived classes must override this function in custom scripts. **
        """
        raise NotImplementedError("LandlabTemplate.initialize_landlab_components() must be overridden in your custom script. " \
                                  "See the doc string in the template file for information on the recommended " \
                                  "way to write this function.")



    def update_until(self, ASPECT_solution_at_Landlab_nodes_dict, ASPECT_additional_info_dict):
        """
        Run the Landlab model for the duration of the ASPECT timestep and return the change in topography at each node.
        This function is where the LandLab components will be called. The change in topography is used to determine the 
        surface velocity of the ASPECT mesh.
        ** Derived classes must override this function in custom scripts. **

        Parameters:
        - ASPECT_solution_at_Landlab_nodes_dict: a dictionary mapping ASPECT variable names to their values at each node on
                                                 the LandLab mesh. This dictionary will include entries for "x velocity" and
                                                 "y velocity" (and "z velocity" for 3D models), as well as the compositional
                                                 fields, and the pressure and temperature.

        - ASPECT_additional_info_dict: a dictionary containing additional information from ASPECT, such as the ASPECT time, 
                                       timestep, dimension, and output directory.
        """

        raise NotImplementedError("LandlabTemplate.update_until() must be overridden in your custom script. " \
                                  "See the doc string in the template file for information on the recommended " \
                                  "way to write this function.")



    # ---------------------------------------------------------------------------
    # The rest of these functions likely will not require any modification.
    # ---------------------------------------------------------------------------
    def determine_uplift_velocity(self, ASPECT_dim, ASPECT_fields_at_Landlab_nodes_dict):
        """
        Determine uplift velocity of the Landlab mesh using the ASPECT velocity. In 3D, the vertical velocity is directly obtained 
        from the z-velocity calculated in ASPECT. In 2D, the vertical velocity is obtained by projecting the y-velocity from the 
        ASPECT surface (which is expected to be located at y=0 on the Landlab mesh) to all nodes on the Landlab mesh.

        For a 2D ASPECT model, the coupling currently only supports structured Landlab grids. The velocity from ASPECT is projected
        across the Landlab grid such that the values from ASPECT are fixed along a given x-value on the Landlab grid.

        Parameters:
        - ASPECT_dim: the dimension of the ASPECT model (2 or 3).
        - ASPECT_fields_at_Landlab_nodes_dict: a dictionary mapping ASPECT variables to values at each node on the Landlab mesh.
        """
        if ASPECT_dim == 2:
            slice_y_velocity = ASPECT_fields_at_Landlab_nodes_dict["y velocity"]

            vertical_velocity = np.zeros(self.model_grid.number_of_nodes)
            unique_x_values = np.unique(self.model_grid.x_of_node)
            for x in unique_x_values:
                vertical_velocity[self.model_grid.x_of_node == x] = slice_y_velocity[unique_x_values == x]

        elif ASPECT_dim == 3:
            vertical_velocity = ASPECT_fields_at_Landlab_nodes_dict["z velocity"]

        return vertical_velocity
    
    def determine_horizontal_velocity(self, ASPECT_dim, ASPECT_fields_at_Landlab_nodes_dict):
        """
        Determine horizontal velocity from ASPECT variables. In 3D, the horizontal velocity is obtained by directly projecting the 
        x and y velocity from ASPECT to the links of the Landlab mesh. In 2D, the horizontal velocity is obtained by projecting the 
        x velocity from ASPECT to the links of the Landlab mesh and setting the y velocity to zero.

        For a 2D ASPECT model, the coupling currently only supports structured Landlab grids. The velocity from ASPECT is projected
        across the Landlab grid such that the values from ASPECT are fixed along a given x-value on the Landlab grid.

        Parameters:
        - ASPECT_dim: the dimension of the ASPECT model (2 or 3).
        - ASPECT_fields_at_Landlab_nodes_dict: a dictionary mapping ASPECT variables to values at each node on the Landlab mesh.
        """
        x_velocity = ASPECT_fields_at_Landlab_nodes_dict["x velocity"]

        if ASPECT_dim == 2:
            projected_x_velocity = np.zeros(self.model_grid.number_of_nodes)
            unique_x_values   = np.unique(self.model_grid.x_of_node)
            for x in unique_x_values:
                projected_x_velocity[self.model_grid.x_of_node == x] = x_velocity[unique_x_values == x]
            
            x_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(projected_x_velocity)
            y_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(np.zeros(self.model_grid.number_of_nodes)) # y velocity is zero since the ASPECT model is 2D.

        elif ASPECT_dim == 3:
            y_velocity = ASPECT_fields_at_Landlab_nodes_dict["y velocity"]

            x_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(x_velocity)
            y_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(y_velocity)

        self.horizontal_velocity[self.model_grid.horizontal_links] = x_vel_at_links[self.model_grid.horizontal_links]
        self.horizontal_velocity[self.model_grid.vertical_links]   = y_vel_at_links[self.model_grid.vertical_links]

        return self.horizontal_velocity

    def dimensional_deposition_erosion(self, ASPECT_dim, deposition_erosion):
        """
        Calculate the change in the topography in a way that is consistent with the dimension expected by the ASPECT model.
        In 3D, this function returns the change in topography at each node. In 2D, this function averages the change in 
        topography across the y-direction and returns the change in topography along y=0, where the ASPECT surface is 
        expected to be located.

        For a 2D ASPECT model, the coupling currently only supports structured Landlab grids. The topography determined by
        Landlab needs to be averaged across the x-direction so that it can be sent to ASPECT as a 1D array.
        
        Parameters:
        - ASPECT_dim: the dimension of the ASPECT model (2 or 3).
        - deposition_erosion: the change in topography at each node on the Landlab mesh
        """
        if ASPECT_dim == 2:
            deposition_erosion_2d = np.zeros(len(np.unique(self.model_grid.x_of_node)))
            unique_x_values = np.unique(self.model_grid.x_of_node)

            for x in unique_x_values:
                deposition_erosion_2d[unique_x_values == x] = np.average(deposition_erosion[self.model_grid.x_of_node == x])
            return deposition_erosion_2d
        
        elif ASPECT_dim == 3:
            return deposition_erosion

    def checkpoint_model_grid(self, checkpoint_dict):
        """
        Checkpoint the Landlab model grid by saving it to a file. This function is called when
        checkpointing the ASPECT model.
        """

        # Extract checkpointing information from the checkpoint dictionary
        checkpoint_index      = checkpoint_dict["Current checkpoint ID"]
        output_directory      = checkpoint_dict["Output directory"]

        # Create LandLab checkpoint directory within the ASPECT output directory.
        output_directory = os.path.join(output_directory, "landlab_checkpoints")
        os.makedirs(output_directory, exist_ok=True)

        print("Checkpointing the LandLab model grid...")
        filename = os.path.join(output_directory, f"landlab_checkpoint_{str(checkpoint_index).zfill(2)}.grid")
        save_grid(self.model_grid, filename, clobber=True)
        pass

    def load_model_grid(self, checkpoint_dict):
        """
        Load the Landlab model grid from a file. This function is called when
        restarting the ASPECT model from a checkpoint.
        """

        restart_checkpoint_id = checkpoint_dict["Resume checkpoint ID"]

        # Extract checkpointing information from the checkpoint dictionary
        output_directory = checkpoint_dict["Output directory"]
        output_directory = os.path.join(output_directory, "landlab_checkpoints")

        # Load the LandLab grid from the checkpoint file corresponding to the checkpoint index.
        print("Loading the LandLab model grid...")
        filename = os.path.join(output_directory, f"landlab_checkpoint_{str(restart_checkpoint_id).zfill(2)}.grid")
        self.model_grid = load_grid(filename)
        self.elevation = self.model_grid.at_node["topographic__elevation"]

        # We need to initialize the components after loading the LandLab grid, since these
        # are not stored.
        self.initialize_landlab_components(None)
        pass

    def get_initial_topography(self, ASPECT_dim):
        """
        Return the initial topography. In 3D, this function returns the initial topography at each node.
        In 2D, this function returns the initial topography along y=0, where the ASPECT surface is expected to be located.

        Parameters:
        - ASPECT_dim: the dimension of the ASPECT model (2 or 3).
        """
        if ASPECT_dim == 2:
            return self.elevation[self.model_grid.y_of_node == 0]
        elif ASPECT_dim == 3:
            return self.elevation
        
    def write_output(self, postprocess_dictionary, output_frequency):
        """
        Write output for visualizing the landlab mesh. This calls a function in the Landlab Python module to write output vtk files. 
        This function is called at the end of each ASPECT timestep after the ASPECT model has been updated and the topography has been evolved.
        The frequency of Landlab output can be controlled by modifying output_frequency.

        Parameters:
        - postprocess_dictionary: a dictionary containing information about the current ASPECT timestep, time, and output directory.
          The dictionary contains the following entries:
            - "ASPECT timestep": the current ASPECT timestep.
            - "ASPECT time": the current ASPECT time.
            - "ASPECT output directory": the directory where ASPECT output is being written.
        """
        time_step  = postprocess_dictionary["ASPECT timestep number"]
        model_time = postprocess_dictionary["ASPECT model time"]
        output_directory = postprocess_dictionary["ASPECT output directory"]
        landlab_output_directory = os.path.join(output_directory, "landlab")

        if time_step % output_frequency != 0:
            return

        if not os.path.isdir(landlab_output_directory):
            os.makedirs(landlab_output_directory, exist_ok=True)

        filename = f"{landlab_output_directory}/landlab_{str(time_step).zfill(5)}.vtk"
        write_legacy_vtk(path=filename, grid=self.model_grid, clobber=True)
        self.vtks.append((model_time, filename))

        with open(f"{landlab_output_directory}/landlab.vtk.series", "w") as f:
            series = {
                "file-series-version": "1.0",
                "files": [
                    {"name": os.path.basename(vtk_name), "time": vtk_time}
                    for vtk_time, vtk_name in self.vtks
                ],
            }
            json.dump(series, f, indent=2)
        pass

    def get_grid_x(self, ASPECT_dim):
        """
        Return the x-coordinates of the grid nodes. In 2D, this function returns the unique x-coordinates.
        In 3D, this function returns the x-coordinates of all nodes. 

        For a 2D ASPECT model, the coupling currently only supports structured Landlab grids. In the case of
        a 2D ASPECT model, the coupling assumes that the ASPECT model is located at y=0 on the Landlab mesh,
        and that the Landlab mesh is structured. Therefore, we only need to extract the unique x-coordinates
        of the Landlab mesh to interface with the ASPECT model.

        
        Parameters:
        - ASPECT_dim: the dimension of the ASPECT model (2 or 3).
        """
        if ASPECT_dim == 2:
            return np.unique(self.model_grid.x_of_node)
        elif ASPECT_dim == 3:
            return self.model_grid.x_of_node
    
    def get_grid_y(self, ASPECT_dim):
        """
        Return the y-coordinates of the grid nodes. In 2D, this function returns an array of zeros equal to 
        the number of unique x-coordinates. In 3D, this function returns the y-coordinates of all nodes.

        For a 2D ASPECT model, the coupling currently only supports structured Landlab grids. In the case of
        a 2D ASPECT model, the coupling assumes that the ASPECT model is located at y=0 on the Landlab mesh,
        and that the Landlab mesh is structured. Therefore, we need to send ASPECT an array of zeros that
        match the length of the number of unique x-coordinates.
        
        Parameters:
        - ASPECT_dim: the dimension of the ASPECT model (2 or 3).
        """
        if ASPECT_dim == 2:
            return np.zeros(np.unique(self.model_grid.x_of_node).size)
        elif ASPECT_dim == 3:
            return self.model_grid.y_of_node
        
    def get_grid_z(self, ASPECT_dim):
        """
        Return the z-coordinates of the grid nodes. This function is only applicable for 3D spherical ASPECT models.

        Parameters:
        - ASPECT_dim: the dimension of the ASPECT model (2 or 3).
        """
        if ASPECT_dim == 3:
            return self.model_grid.z_of_node
        else:
            raise ValueError("get_grid_z is only applicable for 3D ASPECT models.")

    @staticmethod
    def export_ASPECT_functions(subclass, namespace):
        """
        Export the functions of the derived class as module-level functions in the provided namespace. 
        By calling export_ASPECT_functions, the user can simply define a class that inherits LandlabTemplate, override 
        any desired functions, and then call this export function to make the remaining functions available 
        to ASPECT without needing to always define the same boilerplate functions.

        Parameters:
        - subclass: an instance of a class that inherits from LandlabTemplate.
        - namespace: the namespace to which the functions should be added. This is typically done by
                     passing in globals() from the script that defines the derived class.

        Example usage:
        subclass = MyAspectLandlabModel()
        subclass.export_ASPECT_functions(subclass, globals())
        """

        namespace["subclass"] = subclass
        for name in LandlabTemplate.function_explicitly_called_by_aspect:
            namespace[name] = getattr(subclass, name)