/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_deformation_diffusion_h
#define _aspect_mesh_deformation_diffusion_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>


namespace aspect
{
  using namespace dealii;


  namespace MeshDeformation
  {
    /**
     * A plugin that computes the deformation of surface
     * vertices according to the solution of a dim-1 diffusion
     * problem.
     */
    template<int dim>
    class Diffusion : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Diffusion();

        /**
         * Initialize function, which sets the start time and
         * start timestep of diffusion.
         */
        void initialize() override;

        /**
         * The update function sets the current time and determines
         * whether diffusion should be applied in this timestep.
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
         * Declare parameters for the diffusion of the surface.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the diffusion of the surface.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Compute the surface velocity from a difference
         * in surface height given by the solution of
         * the hillslope diffusion problem.
         */
        void diffuse_boundary (const DoFHandler<dim> &free_surface_dof_handler,
                               const IndexSet &mesh_locally_owned,
                               const IndexSet &mesh_locally_relevant,
                               LinearAlgebra::Vector &output,
                               const std::set<types::boundary_id> &boundary_id) const;

        /**
         * Check that the size of the next time step is not larger than the conduction
         * timestep based on the mesh size and the diffusivity on each cell along the
         * surface. The computed time step has to satisfy the CFL
         * number chosen in the input parameter file on each cell of the mesh.
         */
        void check_diffusion_time_step (const DoFHandler<dim> &mesh_deformation_dof_handler,
                                        const std::set<types::boundary_id> &boundary_ids) const;

        /**
         * The hillslope transport coefficient or diffusivity [m2/s]
         * used in the hillslope diffusion of the deformed
         * surface.
         */
        double diffusivity;

        /**
         * Maximum number of steps between the application of diffusion.
         */
        unsigned int timesteps_between_diffusion;

        /**
         * Whether or not in the current timestep, diffusion
         * should be applied.
         */
        bool apply_diffusion;

        /**
         * Boundaries along which the mesh is allowed to move tangentially
         * despite of the Stokes velocity boundary conditions.
         */
        std::set<types::boundary_id> additional_tangential_mesh_boundary_indicators;
    };
  }
}


#endif
