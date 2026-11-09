// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
                                                  "Physical units: $\\text{Pa}$.")
    }
  }
}
