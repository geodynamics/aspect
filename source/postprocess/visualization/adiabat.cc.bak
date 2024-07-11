/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/adiabat.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Adiabat<dim>::
      Adiabat ()
        :
        DataPostprocessor<dim> (),
        Interface<dim>("K,Pa,kg/m/m/m,kg/m/m/m/m")
      {}



      template <int dim>
      std::vector<std::string>
      Adiabat<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;
        solution_names.emplace_back("adiabatic_temperature");
        solution_names.emplace_back("adiabatic_pressure");
        solution_names.emplace_back("adiabatic_density");
        solution_names.emplace_back("adiabatic_density_derivative");
        return solution_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      Adiabat<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(4,
            DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }


      template <int dim>
      UpdateFlags
      Adiabat<dim>::
      get_needed_update_flags () const
      {
        return update_quadrature_points;
      }


      template <int dim>
      void
      Adiabat<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 4,                   ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            computed_quantities[q](0) = this->get_adiabatic_conditions().temperature(input_data.evaluation_points[q]);
            computed_quantities[q](1) = this->get_adiabatic_conditions().pressure(input_data.evaluation_points[q]);
            computed_quantities[q](2) = this->get_adiabatic_conditions().density(input_data.evaluation_points[q]);
            computed_quantities[q](3) = this->get_adiabatic_conditions().density_derivative(input_data.evaluation_points[q]);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Adiabat,
                                                  "adiabat",
                                                  "A visualization output "
                                                  "object that generates "
                                                  "adiabatic temperature, pressure, "
                                                  "density, and density derivative (with regard to depth)"
                                                  "as produced by the \\texttt{AdiabaticConditions} class."
                                                  "\n\n"
                                                  "Physical units: \\si{\\kelvin}, \\si{\\pascal}, "
                                                  "\\si{\\kilo\\gram\\per\\meter\\cubed\\per\\meter}, "
                                                  "respectively, for the four components.")
    }
  }
}
