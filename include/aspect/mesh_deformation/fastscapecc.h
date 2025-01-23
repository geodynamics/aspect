/*
  Copyright (C) 2023 by the authors of the ASPECT code.
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

#ifndef _aspect_mesh_deformation_fastscapecc_h
#define _aspect_mesh_deformation_fastscapecc_h

#include <type_traits>

#include <aspect/global.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/mesh_deformation/fastscapecc_adapter.h>
#include <aspect/simulator_access.h>

#include <fastscapelib/flow/flow_graph.hpp>
#include <fastscapelib/eroders/spl.hpp>

#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>

namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {

    /**
     * A plugin that utilizes the landscape evolution code Fastscapelib (C++)
     * to deform the ASPECT boundary through erosion.
     *
     */
    template<int dim>
    class FastScapecc : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        FastScapecc();

        /**
         * Initialize variables for FastScape.
         */
        void initialize () override;

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
         * Surface mesh and solution.
         */
        using SurfaceMeshType = Triangulation<dim-1, dim>;

        SurfaceMeshType surface_mesh;
        DoFHandler<SurfaceMeshType::dimension, SurfaceMeshType::space_dimension> surface_mesh_dof_handler;
        LinearAlgebra::Vector surface_solution;
        mutable AffineConstraints<double> surface_constraints; // Constraints for hanging nodes
        dealii::LinearAlgebra::distributed::Vector<double> boundary_solution;

        template <class M>
        void init_surface_mesh(M &geom_model);

        template <class M, typename std::enable_if_t<std::is_same_v<M, GeometryModel::Box<3>>>>
        void init_surface_mesh(M &geom_model);

        template <class M, typename std::enable_if_t<std::is_same_v<M, GeometryModel::SphericalShell<3>>>>
        void init_surface_mesh(M &geom_model);

        /**
         * Pointers to Fastscapelib objects
         */
        using GridAdapterType = typename fastscapelib::dealii_grid<SurfaceMeshType>;
        using FlowGraphType = typename fastscapelib::flow_graph<GridAdapterType>;
        std::unique_ptr<GridAdapterType> grid;
        std::unique_ptr<FlowGraphType> flow_graph;
        std::unique_ptr<fastscapelib::spl_eroder<FlowGraphType>> spl_eroder;

        void project_surface_solution(const std::set<types::boundary_id> &boundary_ids);

        /**
         * Project the Stokes velocity solution onto the
         * free surface. Called by make_constraints()
         */
        // void project_velocity_onto_boundary (const DoFHandler<dim> &free_surface_dof_handler,
        //                                      const IndexSet &mesh_locally_owned,
        //                                      const IndexSet &mesh_locally_relevant,
        //                                      LinearAlgebra::Vector &output) const;


        /**
         * Fill velocity data table to be interpolated back onto the ASPECT mesh.
         */
        Table<dim,double> fill_data_table(std::vector<double> &values,
                                          TableIndices<dim> &size_idx,
                                          const int &array_size) const;

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
         * Maximum timestep allowed for FastScape, if the suggested timestep exceeds this
         * limit it is repeatedly divided by 2 until the final timestep is smaller than this parameter.
         */

        /**
         * Check whether FastScape needs to be restarted.
         */
        mutable bool restart;

        /**
         * ASPECT end time to check if we should destroy FastScape.
         */
        double end_time;

        /**
         * Total number of Fastscapelib grid nodes.
         */
        unsigned int n_grid_nodes;

        /**
         * How many levels FastScape should be refined above the maximum ASPECT surface resolution.
         */
        unsigned int additional_refinement_levels;

        /**
         * Maximum expected refinement level at ASPECT's surface.
         * This and resolution_difference are required to properly transfer node data from
         * ASPECT to FastScape.
         */
        int maximum_surface_refinement_level;

        /**
         * Difference in refinement levels expected at the ASPECT surface,
         * where this would be set to 2 if 3 refinement leves are set at the surface.
         * This and surface_resolution are required to properly transfer node data from
         * ASPECT to FastScape.
         *
         * TODO: Should this be kept this way, or make it so the input is the expected levels
         * of refinement at the surface, and we can subtract one within the code? Also,
         * it would be good to find a way to check these are correct, because they are a
         * common source of errors.
         */
        int surface_refinement_difference;

        /**
         * Seed number for initial topography noise in FastScape.
         */
        int fs_seed;

        /**
         * Table for interpolating FastScape surface velocities back to ASPECT.
         */
        std::array< unsigned int, dim > table_intervals;

        /**
         * Magnitude (m) of the initial noise applied to FastScape.
         * Applied as either a + or - value to the topography
         * such that the total difference can be up to 2*noise_h.
         */
        double noise_h;

        // FastScape erosional parameters //
        /**
         * Drainage area exponent for the stream power law (m variable in FastScape surface equation).
         */
        double m;

        /**
         * Slope exponent for the steam power law (n variable in FastScape surface equation).
         */
        double n;

        /**
         * Slope exponent for multi-direction flow, where 0 is uniform, and 10 is steepest descent. (-1 varies with slope)
         * (p variable in FastScape surface equation).
         */
        double p;

        /**
         * Bedrock river incision rate for the stream power law
         * (meters^(1-2m)/yr, kf variable in FastScape surface equation).
         */
        double kff;

        /**
         * Sediment river incision rate for the stream power law (meters^(1-2m)/yr).
         * When set to -1 this is identical to the bedrock value.
         * (kf variable in FastScape surface equation applied to sediment).
         */
        double kfsed;

        /**
         * Bedrock transport coefficient for hillslope diffusion (m^2/yr, kd in FastScape surface equation).
         */
        double kdd;

        /**
         * Bedrock transport coefficient for hillslope diffusion (m^2/yr). When set to -1 this is
         * identical to the bedrock value.
         * (kd in FastScape surface equation applied to sediment).
         */
        double kdsed;

        /**
         * Precision value for how close a ASPECT node must be to the FastScape node
         * for the value to be transferred. This is only necessary if use_v is set to 0
         * and the free surface is used to advect the surface with a normal projection.
         */
        double node_tolerance;

        /**
         * FastScape X extent (ASPECT X extent + 2*dx for ghost nodes).
         */
        double precision;
    };
  }
}

#endif
