/*
  Copyright (C) 2018 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_deformation_boundary_function_h
#define _aspect_mesh_deformation_boundary_function_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {
    template<int dim>
    class BoundaryFunction : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        BoundaryFunction();

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
         * A function object representing the mesh deformation.
         */
        Functions::ParsedFunction<dim> function;
    };
  }
}


#endif
