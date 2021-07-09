/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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

#include <aspect/heating_model/new_heating.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/utilities.h>
#include <aspect/structured_data.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace HeatingModel
  {

    template <int dim>
    NewHeating<dim>::NewHeating ()
      :
      surface_boundary_id(numbers::invalid_unsigned_int)
    {}

    template <int dim>
    void
    NewHeating<dim>::initialize ()
    {
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");
      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(surface_boundary_id);
      Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,
                                                    1);
    }

    template <int dim>
    void
    NewHeating<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                                HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
    
        {
          const double depth = this->get_geometry_model().depth(material_model_inputs.position[q]);
          const double temp = material_model_inputs.temperature[q];
          const double T_surface = this->get_boundary_temperature_manager().minimal_temperature(this->get_fixed_temperature_boundary_indicators());
          const double T_bottom = this->get_boundary_temperature_manager().maximal_temperature(this->get_fixed_temperature_boundary_indicators());
          const double zlo = 125e3;
          double T_profile = 1623;
          //MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
          //MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());
          const double therm_diff = material_model_outputs.thermal_conductivities[q] / (material_model_outputs.densities[q] * material_model_outputs.specific_heat[q]);
          const double heaty_cappy = material_model_outputs.specific_heat[q];
          if (depth > zlo)
            {
              T_profile = T_bottom;
            }
          else
            {
              const double pi = 3.1415;
              const double sea_age = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id,
                                                                                           material_model_inputs.position[q],
                                                                                           0);

              const double exp_fac = -therm_diff * std::pow(pi, 2) * sea_age / std::pow(zlo, 2);
              const double sin_fac = pi / zlo;
              unsigned int n = 1;
              unsigned int m = 5;
              double sum_terms = 0;
              while (n <= m)
                {
                  sum_terms += 1/(double)n * std::exp(std::pow((double)n, 2) * exp_fac) * std::sin((double)n * depth * sin_fac);
                  n += 1;
                }
              T_profile = T_surface + (T_bottom - T_surface) * (depth / zlo + 2 / pi * sum_terms);
            }
          double T_diff = 0;
          if ( (T_profile - temp) > 0 )
            {
              T_diff = T_profile - temp;
            }
          else
            {
              T_diff = 0;
            }
          double dt = this->get_timestep();
          if ( dt == 0)
            {
              dt = 1;
            }
//          heating_model_outputs.heating_source_terms[q] = T_diff * material_model_outputs.densities[q] * heat_cap / dt / 60 / 60 / 24 / 365.25;
//          heating_model_outputs.heating_source_terms[q] = T_diff * material_model_outputs.densities[q] * heaty_cappy / dt / 60 / 60 / 24 / 365.25;
          heating_model_outputs.heating_source_terms[q] = T_diff * material_model_outputs.densities[q] * heaty_cappy / dt;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }

    template <int dim>
    void
    NewHeating<dim>::update ()
    {
      Interface<dim>::update ();
      Utilities::AsciiDataBoundary<dim>::update();
    }

    template <int dim>
    void
    NewHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
//                                                              "$ASPECT_SOURCE_DIR/data/boundary-temperature/ascii-data/test/",
                                                              "/",
                                                              "box_2d_%s.%d.txt");
        prm.enter_subsection("New heating");
        {
          prm.declare_entry ("Heat capacity", "1000",
                             Patterns::Anything(),
                             "The heat capacity in J/K");
          prm.declare_entry ("Thermal diffusivity", "0.8e-7",
                             Patterns::Anything(),
                             "The thermal diffusivity");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    NewHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        Utilities::AsciiDataBoundary<dim>::parse_parameters(prm);

        prm.enter_subsection("New heating");
        {
          heat_cap = prm.get_double("Heat capacity");
          kappa = prm.get_double("Thermal diffusivity");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(NewHeating,
                                  "New heating",
                                  "Implementation of a standard and a simplified model of "
                                  "adiabatic heating.")
  }
}
