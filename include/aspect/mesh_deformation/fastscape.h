/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
 */

#ifndef _aspect_mesh_deformation_fastscape_h
#define _aspect_mesh_deformation_fastscape_h

#include <aspect/global.h>

#ifdef ASPECT_WITH_FASTSCAPE

#include <aspect/mesh_deformation/interface.h>

namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * A plugin that utilizes the landscape evolution code FastScape
     * to deform the ASPECT boundary through advection, uplift,
     * hillslope diffusion, sediment deposition, marine diffusion,
     * and the stream power law, which describes river incision.
     *
     * @ingroup MeshDeformation
     */
    template <int dim>
    class FastScape : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize variables for FastScape.
         */
        virtual void initialize () override;

        /**
         * Destructor for FastScape.
         */
        ~FastScape() override;

        /**
         * A function that creates constraints for the velocity of certain mesh
         * vertices (e.g. the surface vertices) for a specific boundary.
         * The calling class will respect
         * these constraints when computing the new vertex positions.
         */
        virtual
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const override;

        /**
         * Returns whether or not the plugin requires surface stabilization
         */
        bool needs_surface_stabilization () const override;

        /**
         * Declare parameters for the FastScape plugin.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the FastScape plugin.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Function used to set the FastScape ghost nodes. FastScape boundaries are
         * not uplifted or periodic for advection and diffusion. By using a layer
         * of extra nodes in the FastScape model, we can avoid seeing FastScape boundary
         * effects within ASPECT. Similarly, we can use these nodes to have fully
         * periodic boundaries, where we check the flow direction and update the FastScape
         * ghost nodes, and the nodes one layer inward (boundary nodes in ASPECT) to match
         * the parameters on the other side (vx, vy, vz, h). This is done every ASPECT timestep
         * before running FastScape.
         */
        void set_ghost_nodes(std::vector<double> &elevation,
                             std::vector<double> &velocity_x,
                             std::vector<double> &velocity_y,
                             std::vector<double> &velocity_z,
                             const double &fastscape_timestep_in_years,
                             const bool init) const;

        /**
         * Function to determine whether the current index is a ghost node
         */
        bool is_ghost_node(const unsigned int &index,
                           const bool &exclude_boundaries) const;

        /**
         * Function to fill the Fastscape arrays (height and velocities) with the data received from ASPECT in the correct index order.
         */
        void fill_fastscape_arrays(std::vector<double> &elevation,
                                   std::vector<double> &bedrock_transport_coefficient_array,
                                   std::vector<double> &bedrock_river_incision_rate_array,
                                   std::vector<double> &velocity_x,
                                   std::vector<double> &velocity_y,
                                   std::vector<double> &velocity_z,
                                   std::vector<std::vector<double>> &temporary_variables) const;

        /**
         * Function to get the ASPECT topography and velocities at the surface, and an index for transferring these to FastScape.
         */
        std::vector<std::vector<double>> get_aspect_values() const;

        /**
         * Function to initialize or restart FastScape
         */
        void initialize_fastscape(std::vector<double> &elevation,
                                  std::vector<double> &basement,
                                  std::vector<double> &bedrock_transport_coefficient_array,
                                  std::vector<double> &bedrock_river_incision_rate_array,
                                  std::vector<double> &silt_fraction) const;

        /**
         * Execute FastScape
         */
        void execute_fastscape(std::vector<double> &elevation,
                               std::vector<double> &extra_vtk_field,
                               std::vector<double> &velocity_x,
                               std::vector<double> &velocity_y,
                               std::vector<double> &velocity_z,
                               const double &fastscape_timestep_in_years,
                               const unsigned int &fastscape_iterations) const;

        /**
         * Function to apply orographic (mountain related, e.g., wind or elevation)
         * controls to the FastScape model.
         */
        void apply_orographic_controls(const std::vector<double> &elevation,
                                       std::vector<double> &bedrock_river_incision_rate_array,
                                       std::vector<double> &bedrock_transport_coefficient_array) const;

        /**
         * Fill velocity data table to be interpolated back onto the ASPECT mesh.
         */
        Table<dim,double> fill_data_table(std::vector<double> &values,
                                          TableIndices<dim> &size_idx,
                                          const unsigned int &fastscape_nx,
                                          const unsigned int &fastscape_ny) const;

        /**
         * Read data from file for restarting.
         */
        void read_restart_files(std::vector<double> &elevation,
                                std::vector<double> &basement,
                                std::vector<double> &silt_fraction) const;

        /**
         * Save data to file for restarting.
         */
        void save_restart_files(const std::vector<double> &elevation,
                                std::vector<double> &basement,
                                std::vector<double> &silt_fraction) const;

        /**
         * Suggestion for the number of FastScape steps to run for every ASPECT timestep,
         * where the FastScape timestep is determined by ASPECT_timestep_length divided by
         * this parameter.
         */
        unsigned int fastscape_steps_per_aspect_step;

        /**
         * Maximum timestep allowed for FastScape, if the suggested timestep exceeds this
         * limit it is repeatedly divided by 2 until the final timestep is smaller than this parameter.
         */
        double maximum_fastscape_timestep;

        /**
         * Check whether FastScape needs to be restarted. This is used as
         * a mutable bool because we determine whether the model is being resumed in
         * initialize(), and then after reinitializing FastScape we change it to false
         * so it does not initialize FastScape again in future timesteps.
         * TODO: There is probably a better way to do this, and restarts should be rolled into
         * the general ASPECT restart.
         */
        mutable bool restart;

        /**
         * FastScape cell size in X.
         */
        double fastscape_dx;

        /**
         * FastScape cell size in Y.
         */
        double fastscape_dy;

        /**
         * FastScape X extent (ASPECT X extent + 2*dx for ghost nodes).
         */
        double fastscape_x_extent;

        /**
         * Fastscape Y extent (ASPECT Y extent + 2*dy for ghost nodes).
         */
        double fastscape_y_extent;

        /**
         * User set FastScape Y extent for a 2D ASPECT model.
         */
        double fastscape_y_extent_2d;

        /**
         * Number of x points in FastScape array.
         */
        unsigned int fastscape_nx;

        /**
         * Number of y points in FastScape array.
         */
        unsigned int fastscape_ny;

        /**
         * Vertical exaggeration in FastScape visualization.
         */
        double vexp;

        /**
         * How many levels FastScape should be refined above the maximum ASPECT surface resolution.
         */
        unsigned int additional_refinement_levels;

        /**
         * Maximum expected refinement level at ASPECT's surface.
         * This and resolution_difference are required to properly transfer node data from
         * ASPECT to FastScape.
         */
        unsigned int maximum_surface_refinement_level;

        /**
         * Difference in refinement levels expected at the ASPECT surface,
         * where this would be set to 2 if 3 refinement levels are set at the surface.
         * This and surface_resolution are required to properly transfer node data from
         * ASPECT to FastScape.
         *
         * TODO: Should this be kept this way, or make it so the input is the expected levels
         * of refinement at the surface, and we can subtract one within the code? Also,
         * it would be good to find a way to check these are correct, because they are a
         * common source of errors.
         */
        unsigned int surface_refinement_difference;

        /**
         * If set to true, the FastScape surface is averaged along Y and returned
         * to ASPECT. If set to false, the center slice of the FastScape model is
         * returned to ASPECT.
         */
        bool average_out_of_plane_surface_topography;

        /**
         * Seed number for initial topography noise in FastScape.
         */
        int fastscape_seed;

        /**
         * Variable to hold ASPECT domain extents.
         */
        std::array<std::pair<double,double>,dim> grid_extent;

        /**
         * Table for interpolating FastScape surface velocities back to ASPECT.
         */
        std::array<unsigned int, dim> table_intervals;

        /**
         * Whether or not to use the ghost nodes.
         */
        bool use_ghost_nodes;

        /**
         * Magnitude (m) of the initial noise applied to FastScape.
         * Applied as either a + or - value to the topography
         * such that the total difference can be up to 2*noise_elevation.
         */
        double noise_elevation;

        /**
         * Sediment rain in m/yr, added as a flat increase to the FastScape surface
         * in the marine domain every ASPECT timestep before running FastScape.
         */
        std::vector<double> sediment_rain_rates;

        /**
         * Time at which each interval of sediment_rain_rates is active. Should contain
         * one less value than sediment_rain_rates, assuming that sediment_rain_rates[0]
         * is applied from model time 0 until sediment_rain_times[0].
         */
        std::vector<double> sediment_rain_times;

        /**
         * Flag for having FastScape advect/uplift the surface. If the free surface is used
         * in conjunction with FastScape, this can be set to false, then FastScape will only
         * apply erosion/deposition to the surface and not advect or uplift it.
         */
        bool fastscape_advection_uplift;

        /**
         * Node tolerance for how close a ASPECT node must be to the FastScape node
         * for the value to be transferred. This is only necessary if use_v is set to 0
         * and the free surface is used to advect the surface with a normal projection, or
         * if there is a surface refinement level difference leading to excess interpolation
         * points in areas of high ASPECT resolution.
         */
        double node_tolerance;

        /**
         * Interval between the generation of graphical output. This parameter
         * is read from the input file and consequently is not part of the
         * state that needs to be saved and restored.
         */
        double output_interval;

        /**
         * A time (in seconds) at which the last graphical output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        mutable double last_output_time;

        /**
         * @name Fastscape boundary conditions
         * @{
         */

        /**
         * FastScape bottom boundary condition that determines topography at the FastScape bottom boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., top and bottom are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int bottom;

        /**
         * FastScape top boundary condition that determines topography at the FastScape top boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., top and bottom are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int top;

        /**
         * FastScape right boundary condition that determines topography at the FastScape right boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., left and right are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int right;

        /**
         * FastScape left boundary condition that determines topography at the FastScape left boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., left and right are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int left;

        /**
         * Parameters that set the FastScape boundaries periodic even though the ghost nodes are set 'fixed'
         */
        bool topbottom_ghost_nodes_periodic;
        bool leftright_ghost_nodes_periodic;

        /**
         * Integer that holds the full boundary conditions sent to FastScape (e.g., 1111).
         */
        unsigned int fastscape_boundary_conditions;

        /**
         * Prescribed flux per unit length into the model through the bottom boundary (m^2/yr).
         */
        double bottom_flux;

        /**
         * Prescribed flux per unit length into the model through the top boundary (m^2/yr).
         */
        double top_flux;

        /**
         * Prescribed flux per unit length into the model through the right boundary (m^2/yr).
         */
        double right_flux;

        /**
         * Prescribed flux per unit length into the model through the left boundary (m^2/yr).
         */
        double left_flux;
        /**
         * @}
         */

        /**
         * @name FastScape subaerial erosional parameters
         * @{
         */

        /**
         * Drainage area exponent for the stream power law. ($m$ variable in FastScape surface equation.)
         */
        double drainage_area_exponent_m;

        /**
         * Slope exponent for the stream power law. ($n$ variable in FastScape surface equation.)
         */
        double slope_exponent_n;

        /**
         * Slope exponent for multi-direction flow, where 0 is uniform, and 10 is steepest descent. (-1 varies with slope.)
         * ($p$ variable in FastScape surface equation.)
         */
        double slope_exponent_p;

        /**
         * Bedrock deposition coefficient. Higher values deposit more sediment
         * inside the domain. ($G$ variable in FastScape surface equation.)
         */
        double bedrock_deposition_g;

        /**
         * Sediment deposition coefficient. Higher values deposit more sediment inside the domain.
         * When set to -1 this is identical to the bedrock value.
         * ($G$ variable in FastScape surface equation applied to sediment.)
         */
        double sediment_deposition_g;

        /**
         * Bedrock river incision rate for the stream power law.
         * (meters^(1-2m)/yr, $kf$ variable in FastScape surface equation.)
         */
        double bedrock_river_incision_rate;

        /**
         * Sediment river incision rate for the stream power law (meters^(1-2m)/yr).
         * When set to -1 this is identical to the bedrock value.
         * ($kf$ variable in FastScape surface equation applied to sediment.)
         */
        double sediment_river_incision_rate;

        /**
         * Bedrock transport coefficient for hillslope diffusion (m^2/yr, kd in FastScape surface equation.)
         */
        double bedrock_transport_coefficient;

        /**
         * Sediment transport coefficient for hillslope diffusion (m^2/yr). When set to -1 this is
         * identical to the bedrock value.
         * (kd in FastScape surface equation applied to sediment).
         */
        double sediment_transport_coefficient;
        /**
         * @}
         */

        /**
         * @name Fastscape marine parameters
         * @{
         */

        /**
         * Fastscape sea level (m), set relative to the ASPECT surface where
         * a sea level of zero will represent the initial maximum unperturbed
         * Y (2D) or Z (3D) extent of the ASPECT domain. A negative value of
         * the sea level means the sea level lies below the initial unperturbed
         * top boundary of the domain.
         */
        double sea_level;

        /**
         * Parameters to set an extra erosional base level
         * on the ghost nodes that differs from sea level.
         */
        bool use_fixed_erosional_base;

        /**
         * Height of the extra erosional base level.
         */
        double h_erosional_base;

        /**
         * Surface porosity for sand.
         */
        double sand_surface_porosity;

        /**
         * Surface porosity for silt.
         */
        double silt_surface_porosity;

        /**
         * Sands e-folding depth for exponential porosity law (m).
         */
        double sand_efold_depth;

        /**
         * Silts e-folding depth for exponential porosity law (m).
         */
        double silt_efold_depth;

        /**
         * Sand-silt ratio
         */
        double sand_silt_ratio;

        /**
         * Averaging depth/thickness for sand-silt equation (m).
         */
        double sand_silt_averaging_depth;

        /**
         * Sand marine transport coefficient. (marine diffusion, m^2/yr.)
         */
        double sand_transport_coefficient;

        /**
         * Silt marine transport coefficient. (marine diffusion, m^2/yr.)
         */
        double silt_transport_coefficient;

        /**
         * Flag to use the marine component of FastScape.
         */
        bool use_marine_component;
        /**
         * @}
         */

        /**
         * @name Orographic parameters
         * @{
         */

        /**
         * Set a flat height (m) after which the flat_erosional_factor
         * is applied to the bedrock river incision rate and transport coefficient.
         */
        int flat_elevation;

        /**
         * Set the height (m) after which the model will track the ridge line
         * and based on the wind direction will apply the wind_barrier_erosional_factor
         * to the bedrock river incision rate and transport coefficient.
         */
        int wind_barrier_elevation;

        /**
         * Wind direction for wind_barrier_erosional_factor.
         */
        unsigned int wind_direction;

        /**
         * Factor to multiply the bedrock river incision rate and transport coefficient by depending on them
         * flat_elevation.
         */
        double flat_erosional_factor;

        /**
         * Factor to multiply the bedrock river incision rate and transport coefficient by depending on them
         * wind_barrier_elevation and wind direction.
         */
        double wind_barrier_erosional_factor;

        /**
         * Flag to stack both orographic controls.
         */
        bool stack_controls;

        /**
         * Flag to use orographic controls.
         */
        bool use_orographic_controls;
        /**
         * @}
         */
    };
  }
}

#endif
#endif
