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


#ifndef _aspect_postprocess_surface_topography_h
#define _aspect_postprocess_surface_topography_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes dynamic topography at the top and bottom of the domain.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class DynamicTopography : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for the dynamic topography.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Return the topography vector as calculated by the consistent
         * boundary flux (CBF) formulation.
         * The velocity components store the surface-normal traction,
         * and the temperature component stores the dynamic topography
         * computed from that traction.
         */
        const LinearAlgebra::BlockVector &
        topography_vector() const;

        /**
         * Return the cell-wise topography vector as calculated by the CBF formulation,
         * where indices of the vector correspond to cell indices.
         * This vector is considerably smaller than the full topography vector returned
         * by topography_vector(), and is useful for text output and visualization.
         */
        const Vector<float> &
        cellwise_topography() const;

        /**
         * Register the other postprocessor that we need: BoundaryPressures
         */
        std::list<std::string>
        required_other_postprocessors() const override;

        /**
         * Parse the parameters for the postprocessor.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Declare the parameters for the postprocessor.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

      private:
        /**
         * Output the dynamic topography solution to
         * a file for a given boundary id. All the values
         * in @p position_and_topography are written to a file
         * with a file name determined by @p boundary_id.
         */
        void output_to_file(const types::boundary_id boundary_id,
                            const std::vector<std::pair<Point<dim>, double>> &position_and_topography);

        /**
         * A vector which stores the surface stress values calculated
         * by the postprocessor.
         */
        LinearAlgebra::BlockVector topo_vector;

        /**
         * A vector which stores the surface stress values calculated
         * at the midpoint of each surface cell face. This can be
         * given to a visualization postprocessor to output dynamic topography.
         */
        Vector<float> visualization_values;

        /**
         * A parameter that allows users to set the density value
         * above the top surface.
         */
        double density_above;

        /**
         * A parameter that allows users to set the density value
         * below the bottom surface.
         */
        double density_below;

        /**
         * Whether to output the surface topography.
         */
        bool output_surface;

        /**
         * Whether to output the bottom topography.
         */
        bool output_bottom;
    };
  }
}


#endif
