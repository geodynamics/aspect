/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_heat_flux_map_h
#define _aspect_postprocess_visualization_heat_flux_map_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A postprocessor that computes the heat flux density through the top and bottom boundaries.
       *
       * @ingroup Postprocessing
       */
      template <int dim>
      class HeatFluxMap
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          HeatFluxMap();

          /**
           * Initialize the postprocessor.
           */
          void
          initialize() override;

          /**
           * Fill the temporary storage variables with the
           * heat flux for the current time step.
           *
           * @copydoc Interface<dim>::update()
           */
          void update() override;

          /**
           * Compute the heat flux for the given input cell.
           *
           * @copydoc dealii::DataPostprocessor<dim>::evaluate_vector_field()
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * @copydoc Interface<dim>::declare_parameters()
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * @copydoc Interface<dim>::parse_parameters()
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * A flag that determines whether to use the point-wise
           * heat flux calculation or the cell-wise averaged calculation.
           */
          bool output_point_wise_heat_flux;

          /**
           * A temporary storage place for the point-wise heat flux
           * solution. Only initialized and used if output_point_wise_heat_flux
           * is set to true.
           */
          LinearAlgebra::BlockVector heat_flux_density_solution;

          /**
           * A temporary storage place for the cell-wise heat flux
           * solution. Only initialized and used if output_point_wise_heat_flux
           * is set to false.
           */
          std::vector<std::vector<std::pair<double, double>>> heat_flux_and_area;
      };
    }
  }
}


#endif
