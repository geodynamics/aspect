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


#ifndef _aspect_postprocess_heat_flux_map_h
#define _aspect_postprocess_heat_flux_map_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      /**
       * Compute the heat flux for boundaries with prescribed temperature (Dirichlet
       * boundary conditions) using the consistent boundary flux method. The method
       * is described in
       *
       * Gresho, P. M., Lee, R. L., Sani, R. L., Maslanik, M. K., & Eaton, B. E. (1987).
       * The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids
       * problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.
       *
       * The implementation of the method is benchmarked in
       *
       * Juliane Dannberg, Rene Gassmöller, Daniele Thallner, Frederick LaCombe, Courtney Sprain (2024).
       * Changes in core–mantle boundary heat flux patterns throughout the supercontinent cycle,
       * Geophysical Journal International, Volume 237, Issue 3, June 2024, Pages 1251–1274,
       * https://doi.org/10.1093/gji/ggae075.
       *
       * In summary, the method solves the temperature equation again on the boundary faces, with known
       * temperatures and solving for the boundary fluxes that satisfy the equation. Since the
       * equation is only formed on the faces and it can be solved using only diagonal matrices,
       * the computation is cheap. Conceptually simpler methods like evaluating the temperature
       * gradient on the face are significantly less accurate.
       *
       * The function returns a solution vector, which contains the heat flux in the temperature
       * block of the vector.
       */
      template <int dim>
      LinearAlgebra::BlockVector
      compute_dirichlet_boundary_heat_flux_solution_vector (const SimulatorAccess<dim> &simulator_access);

      /**
       * This function computes the combined heat flux through each boundary face (conductive + advective).
       * For reflecting boundaries the conductive heat flux is 0, for boundaries with prescribed heat flux
       * (inhomogeneous Neumann boundary conditions) it is simply the integral of the prescribed heat flux over
       * the face, and for boundaries with non-tangential velocities the advective heat flux is computed as
       * the integral over the advective heat flux density.
       * For boundaries with prescribed temperature (Dirichlet boundary conditions) the heat flux
       * is computed using the compute_dirichlet_boundary_heat_flux_solution_vector() function.
       *
       * The function returns a vector with as many entries as active cells. For each locally owned
       * cell it contains a vector with one entry per face. Each of these entries contains a pair
       * of doubles, containing the combined heat flux (first entry) and face area (second entry).
       * This function is a helper function that unifies the complex heat flux computation necessary
       * for several postprocessors.
       */
      template <int dim>
      std::vector<std::vector<std::pair<double, double>>>
      compute_heat_flux_through_boundary_faces (const SimulatorAccess<dim> &simulator_access);
    }

    /**
     * A postprocessor that computes the point-wise heat flux density through the boundaries.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class HeatFluxMap : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the postprocessor.
         */
        void
        initialize() override;

        /**
         * Evaluate the solution for the heat flux.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

      private:
        /**
         * Output the heat flux density for the boundary determined
         * by @p boundary_id to a file. The heat flux density is
         * handed over in the vector @p heat_flux_and_area. This vector
         * is expected to be of the structure described for the return value
         * of the function compute_heat_flux_through_boundary_faces() and
         * only the values at the faces of the given @p boundary_id are
         * written to the file.
         */
        void output_to_file(const types::boundary_id boundary_id,
                            const std::vector<std::vector<std::pair<double, double>>> &heat_flux_and_area);
    };
  }
}


#endif
