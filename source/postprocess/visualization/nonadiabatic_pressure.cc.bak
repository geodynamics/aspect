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


#include <aspect/postprocess/visualization/nonadiabatic_pressure.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      NonadiabaticPressure<dim>::
      NonadiabaticPressure ()
        :
        DataPostprocessorScalar<dim> ("nonadiabatic_pressure",
                                      update_values | update_quadrature_points),
        Interface<dim>("Pa")
      {}



      template <int dim>
      void
      NonadiabaticPressure<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components, ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double pressure    = input_data.solution_values[q][this->introspection().component_indices.pressure];

            computed_quantities[q](0) = pressure - this->get_adiabatic_conditions().pressure(input_data.evaluation_points[q]);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(NonadiabaticPressure,
                                                  "nonadiabatic pressure",
                                                  "A visualization output object that generates output "
                                                  "for the non-adiabatic component of the pressure."
                                                  "\n\n"
                                                  "The variable that is outputted this way is computed by "
                                                  "taking the pressure at each point and subtracting "
                                                  "from it the adiabatic pressure computed at the beginning "
                                                  "of the simulation. Because the adiabatic pressure is "
                                                  "one way of defining a static pressure background field, "
                                                  "what this visualization postprocessor therefore produces is "
                                                  "\\textit{one} way to compute a \\textit{dynamic "
                                                  "pressure}. There are, however, other ways as well, "
                                                  "depending on the choice of the ``background pressure''."
                                                  "\n\n"
                                                  "Physical units: \\si{\\pascal}.")
    }
  }
}
