/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#include <aspect/postprocess/visualization/mantle_flux.h>

#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/postprocess/mantle_flux_statistics.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      MantleFlux<dim>::MantleFlux ()
        :
        DataPostprocessor<dim> (),
        Interface<dim>("")
      {}



      template <int dim>
      std::vector<std::string>
      MantleFlux<dim>::get_names () const
      {
        return {"mantle_structure",
                "temperature_anomaly_flux_density",
                "thermal_buoyancy_mass_flux_density",
                "thermal_buoyancy_force_rate_density"
               };
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      MantleFlux<dim>::get_data_component_interpretation () const
      {
        return std::vector<DataComponentInterpretation::DataComponentInterpretation>
               (get_names().size(), DataComponentInterpretation::component_is_scalar);
      }



      template <int dim>
      std::string
      MantleFlux<dim>::get_physical_units () const
      {
        return this->convert_output_to_years()
               ?
               "-,K m/year,kg/(m^2 year),N/(m^2 year)"
               :
               "-,K m/s,kg/(m^2 s),N/(m^2 s)";
      }



      template <int dim>
      UpdateFlags
      MantleFlux<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_values | update_quadrature_points;
      }



      template <int dim>
      void
      MantleFlux<dim>::evaluate_vector_field (
        const DataPostprocessorInputs::Vector<dim> &input_data,
        std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,
                ExcInternalError());
        Assert (computed_quantities[0].size() == get_names().size(),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,
                ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> material_inputs(input_data,
                                                                this->introspection());
        MaterialModel::MaterialModelOutputs<dim> material_outputs(n_quadrature_points,
                                                                  this->n_compositional_fields());
        material_inputs.requested_properties =
          MaterialModel::MaterialProperties::density |
          MaterialModel::MaterialProperties::thermal_expansion_coefficient;
        this->get_material_model().evaluate(material_inputs, material_outputs);

        const auto &mantle_flux_statistics =
          this->get_postprocess_manager().template
          get_matching_active_plugin<Postprocess::MantleFluxStatistics<dim>>();

        const double time_scaling =
          this->convert_output_to_years() ? year_in_seconds : 1.0;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const Point<dim> &position = material_inputs.position[q];
            const Tensor<1,dim> gravity =
              this->get_gravity_model().gravity_vector(position);
            const auto values = mantle_flux_statistics.evaluate_point(
                                  position,
                                  material_inputs.velocity[q],
                                  material_inputs.temperature[q],
                                  material_outputs.densities[q],
                                  material_outputs.thermal_expansion_coefficients[q],
                                  gravity);

            computed_quantities[q][0] = values.structure;
            computed_quantities[q][1] = values.temperature_flux_density * time_scaling;
            computed_quantities[q][2] =
              values.buoyancy_mass_flux_density * time_scaling;
            computed_quantities[q][3] =
              values.buoyancy_force_rate_density * time_scaling;
          }
      }



      template <int dim>
      std::list<std::string>
      MantleFlux<dim>::required_other_postprocessors () const
      {
        return {"mantle flux statistics"};
      }
    }
  }
}



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(
        MantleFlux,
        "mantle flux",
        "Writes four fields used by the mantle flux statistics postprocessor: "
        "the detected mantle structure, temperature-anomaly flux density, and "
        "the thermal-buoyancy mass-flux and force-rate densities. Radial velocity and nonadiabatic "
        "temperature are available through existing visualization postprocessors. "
        "The mantle structure field is 1 for plume material, -1 for slab "
        "material, and 0 for background mantle.")
    }
  }
}
