/*
  Copyright (C) 2018 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_deformation_openlem_h
#define _aspect_mesh_deformation_openlem_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>
#include <openlem.cpp>
#include <deal.II/base/parsed_function.h>

#include <aspect/mesh_deformation/openlem.h>
namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * A class that represents a mesh deformation function that can be
     * prescribed on the boundary of the domain.
     *
     * @ingroup MeshDeformation
     */
    template <int dim>
    class OpenLEM : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        OpenLEM();

        /**
         * Initialize variables for openlem.
         */
        virtual void initialize () override;

        /**
         *
         */
        void update() override;

        /**
         * A function that creates constraints for the velocity of certain mesh
         * vertices (e.g. the surface vertices) for a specific boundary.
         * The calling class will respect
         * these constraints when computing the new vertex positions.
         */
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const override;

        /**
         * Returns whether or not the plugin requires surface stabilization
         */
        bool needs_surface_stabilization () const override;

        /**
         * Declare parameters for the free surface handling.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the free surface handling.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Execute openlem
         */
        void execute_openlem(std::vector<double> &elevation,
                               std::vector<double> &extra_vtk_field,
                               std::vector<double> &velocity_x,
                               std::vector<double> &velocity_y,
                               std::vector<double> &velocity_z,
                               const double &openlem_timestep_in_years,
                               const unsigned int &openlem_iterations) const;

        /**
         * Function to fill the openlem arrays (height and velocities) with the data received from ASPECT in the correct index order.
         */
        void fill_openlem_arrays(openlem::Grid<openlem::Node> &grid,
				   // std::vector<double> &elevation,
                                   // std::vector<double> &bedrock_transport_coefficient_array,
                                   // std::vector<double> &bedrock_river_incision_rate_array,
                                   // std::vector<double> &velocity_x,
                                   // std::vector<double> &velocity_y,
                                   // std::vector<double> &velocity_z,
                                   std::vector<std::vector<double>> &temporary_variables);

        /**
         * Function to get the ASPECT topography and velocities at the surface, and an index for transferring these to openlem.
         */
        std::vector<std::vector<double>> get_aspect_values() const;
        
	/**
         * define a old and new grid, so that we can compute a difference
         */
	openlem::Grid<> grid_old;
	openlem::Grid<> grid_new;
	unsigned int openlem_nx;
	unsigned int openlem_ny;
	unsigned int openlem_dx; 
	unsigned int openlem_dy;
	unsigned int openlem_x_extent; 
	unsigned int openlem_y_extent; 
        /**
         * Variable to hold ASPECT domain extents.
         */
        std::array<std::pair<double,double>,dim> grid_extent;

        /**
         * Table for interpolating openlem surface velocities back to ASPECT.
         */
        std::array<unsigned int, dim> table_intervals;

        /**
         * Whether or not to use the ghost nodes.
         */
        bool use_ghost_nodes;

        /**
         * Check whether openlem needs to be restarted. This is used as
         * a mutable bool because we determine whether the model is being resumed in
         * initialize(), and then after reinitializing openlem we change it to false
         * so it does not initialize openlem again in future timesteps.
         * TODO: There is probably a better way to do this, and restarts should be rolled into
         * the general ASPECT restart.
         */
        mutable bool restart;

        /**
         * How many levels openlem should be refined above the maximum ASPECT surface resolution.
         */
        unsigned int additional_refinement_levels;

        /**
         * Maximum expected refinement level at ASPECT's surface.
         * This and resolution_difference are required to properly transfer node data from
         * ASPECT to openlem.
         */
        unsigned int maximum_surface_refinement_level;

        /**
         * User set openlem Y extent for a 2D ASPECT model.
         */
        double openlem_y_extent_2d;

        /**
         * A time (in seconds) at which the last graphical output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        mutable double last_output_time;

        /**
         * Suggestion for the number of openlem steps to run for every ASPECT timestep,
         * where the openlem timestep is determined by ASPECT_timestep_length divided by
         * this parameter.
         */
        unsigned int openlem_steps_per_aspect_step;

        /**
         * Maximum timestep allowed for openlem, if the suggested timestep exceeds this
         * limit it is repeatedly divided by 2 until the final timestep is smaller than this parameter.
         */
        double maximum_openlem_timestep;

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
         * Node tolerance for how close a ASPECT node must be to the FastScape node
         * for the value to be transferred. This is only necessary if use_v is set to 0
         * and the free surface is used to advect the surface with a normal projection, or
         * if there is a surface refinement level difference leading to excess interpolation
         * points in areas of high ASPECT resolution.
         */
        double node_tolerance;
    };
  }
}


#endif
