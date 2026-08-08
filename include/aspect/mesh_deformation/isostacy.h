/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE. If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_mesh_deformation_isostacy_h
#define _aspect_mesh_deformation_isostacy_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * Compute an initial mesh deformation that balances the mass of
     * vertical columns and produces an isostatically compensated
     * surface.
     */
    template <int dim>
    class Isostacy
      :
      public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:

        Isostacy ();

        /**
         * Compute the initial isostatic deformation.
         */
        void
        initialize () override;

        /**
         * Return the prescribed deformation on a boundary.
         */
        Tensor<1,dim>
        compute_initial_deformation_on_boundary (
          const types::boundary_id boundary_indicator,
          const Point<dim> &position) const override;

        /**
         * Whether the mesh deformation requires free surface stabilization.
         */
        bool
        needs_surface_stabilization () const override;

        /**
         * Declare parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * The number of equally spaced lateral sampling points used to
         * sample the mesh surface. For each sampling point, a vertical
         * column is constructed through the model domain to compute the
         * column mass.
         */
        unsigned int n_lateral_points;

        /**
         * The number of equally spaced vertical sampling points used
         * within each column. These points are used to sample the
         * material properties along the column and numerically
         * integrate the column mass.
         */
        unsigned int n_vertical_points;

        /**
         * The depth above which the model is assumed to be in
         * sostatic equilibrium.
         */
        double compensation_depth;

        /**
         * The characteristic horizontal length scale (in m) over which
         * topographic loads are assumed to attain isostatic compensation.
         */
        double isostatic_length_scale;

        /**
         * Maximum allowed initial topography.
         */
        double max_isostatic_topography;

        /**
         * Computed topography profile.
         */
        std::vector<double> topography;
    };
  }
}

#endif
