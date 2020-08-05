/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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



#include <aspect/heating_model/radioactive_decay.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/global.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    RadioactiveDecay<dim>::RadioactiveDecay ()
    {}


    template <int dim>
    void
    RadioactiveDecay<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      AssertThrow(crust_composition_num < material_model_inputs.composition[0].size(),
                  ExcMessage("The composition number of crust is larger than number of composition fields."));

      for (unsigned int q = 0; q < heating_model_outputs.heating_source_terms.size(); ++q)
        {
          double timedependent_radioactive_heating_rates = 0;

          if (n_radio_heating_elements != 0)
            {
              double crust_fraction = 0;

              if (is_crust_defined_by_composition)
                {
                  crust_fraction = material_model_inputs.composition[q][crust_composition_num];

                  if (crust_fraction < 0.0)
                    crust_fraction = 0;
                  if (crust_fraction > 1.0)
                    crust_fraction = 1;
                }
              else if ((this->get_geometry_model()).depth(material_model_inputs.position[q]) < crust_depth)
                crust_fraction = 1;

              for (unsigned element = 0; element < n_radio_heating_elements; ++element)
                {
                  timedependent_radioactive_heating_rates += radioactive_heating_rates[element]
                                                             * (radioactive_initial_concentrations_mantle[element] * (1-crust_fraction)
                                                                + radioactive_initial_concentrations_crust[element] * crust_fraction)
                                                             * std::pow(0.5,this->get_time()/half_decay_times[element]);
                }
            }

          heating_model_outputs.heating_source_terms[q] = timedependent_radioactive_heating_rates
                                                          * material_model_outputs.densities[q];

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }


    template <int dim>
    void
    RadioactiveDecay<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Radioactive decay");
        {
          prm.declare_entry("Number of elements","0",
                            Patterns::Integer(0),
                            "Number of radioactive elements");
          prm.declare_entry("Heating rates","",
                            Patterns::List (Patterns::Double ()),
                            "Heating rates of different elements (W/kg)");
          prm.declare_entry("Half decay times","",
                            Patterns::List (Patterns::Double (0.)),
                            "Half decay times. Units: (Seconds), or "
                            "(Years) if set `use years instead of seconds'.");
          prm.declare_entry("Initial concentrations crust","",
                            Patterns::List (Patterns::Double (0.)),
                            "Initial concentrations of different elements (ppm)");
          prm.declare_entry("Initial concentrations mantle","",
                            Patterns::List (Patterns::Double (0.)),
                            "Initial concentrations of different elements (ppm)");
          prm.declare_entry("Crust defined by composition","false",
                            Patterns::Bool(),
                            "Whether crust defined by composition or depth");
          prm.declare_entry("Crust depth","0.",
                            Patterns::Double(),
                            "Depth of the crust when crust if defined by depth. "
                            "Units: \\si{\\meter}.");
          prm.declare_entry("Crust composition number","0",
                            Patterns::Integer(0),
                            "Which composition field should be treated as crust");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    RadioactiveDecay<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Radioactive decay");
        {

          n_radio_heating_elements  = prm.get_integer ("Number of elements");

          radioactive_heating_rates = Utilities::string_to_double
                                      (Utilities::split_string_list
                                       (prm.get("Heating rates")));

          half_decay_times          = Utilities::string_to_double
                                      (Utilities::split_string_list
                                       (prm.get("Half decay times")));

          radioactive_initial_concentrations_crust  = Utilities::string_to_double
                                                      (Utilities::split_string_list
                                                       (prm.get("Initial concentrations crust")));

          radioactive_initial_concentrations_mantle = Utilities::string_to_double
                                                      (Utilities::split_string_list
                                                       (prm.get("Initial concentrations mantle")));

          is_crust_defined_by_composition = prm.get_bool    ("Crust defined by composition");
          crust_depth                     = prm.get_double  ("Crust depth");
          crust_composition_num           = prm.get_integer ("Crust composition number");


          AssertThrow(radioactive_heating_rates.size() == n_radio_heating_elements,
                      ExcMessage("Number of heating rate entities does not match "
                                 "the number of radioactive elements."));
          AssertThrow(half_decay_times.size() == n_radio_heating_elements,
                      ExcMessage("Number of half decay time entities does not match "
                                 "the number of radioactive elements."));
          AssertThrow(radioactive_initial_concentrations_crust.size() == n_radio_heating_elements,
                      ExcMessage("Number of initial concentration entities in crust "
                                 "does not match the number of radioactive elements."));
          AssertThrow(radioactive_initial_concentrations_mantle.size() == n_radio_heating_elements,
                      ExcMessage("Number of initial concentration entities in mantle "
                                 "does not match the number of radioactive elements."));

          // if we get half_decay_times passed as years convert it to seconds
          if (this->convert_output_to_years())
            for (unsigned int i = 0; i < n_radio_heating_elements; ++i)
              half_decay_times[i] *= year_in_seconds;

          // Convert ppm to SI concentrations
          for (unsigned int i = 0; i < n_radio_heating_elements; ++i)
            {
              radioactive_initial_concentrations_crust[i] *= 1e-6;
              radioactive_initial_concentrations_mantle[i] *= 1e-6;
            }
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
    ASPECT_REGISTER_HEATING_MODEL(RadioactiveDecay,
                                  "radioactive decay",
                                  "Implementation of a model in which the internal "
                                  "heating rate is radioactive decaying in the following rule:\n"
                                  "\\[(\\text{initial concentration})\\cdot 0.5^{\\text{time}/(\\text{half life})}\\]\n"
                                  "The crust and mantle can have different concentrations, and the crust can be "
                                  "defined either by depth or by a certain compositional field.\n"
                                  "The formula is interpreted as having units W/kg.")
  }
}
