/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/nonadiabatic_temperature.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      NonadiabaticTemperature<dim>::
      NonadiabaticTemperature ()
        :
        DataPostprocessorScalar<dim> ("nonadiabatic_temperature",
                                      update_values | update_quadrature_points),
        Interface<dim>("K")
      {}



      template <int dim>
      void
      NonadiabaticTemperature<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double temperature = input_data.solution_values[q][this->introspection().component_indices.temperature];

            computed_quantities[q](0) = temperature - this->get_adiabatic_conditions().temperature(input_data.evaluation_points[q]);
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(NonadiabaticTemperature,
                                                  "nonadiabatic temperature",
                                                  "A visualization output object that generates output "
                                                  "for the non-adiabatic component of the temperature."
                                                  "\n\n"
                                                  "Physical units: \\si{\\kelvin}.")
    }
  }
}
