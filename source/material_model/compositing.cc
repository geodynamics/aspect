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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <deal.II/base/std_cxx11/array.h>
#include <aspect/material_model/compositing.h>
#include <utility>
#include <limits>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    // Parse property names
    template <int dim>
    Property::MaterialProperty
    Compositing<dim>::parse_property_name(const std::string &s)
    {
      if (s == "viscosity")
        return Property::viscosity;
      else if (s == "density")
        return Property::density;
      else if (s == "thermal expansion coefficient")
        return Property::thermal_expansion_coefficient;
      else if (s == "specific heat")
        return Property::specific_heat;
      else if (s == "compressibility")
        return Property::compressibility;
      else if (s == "entropy derivative pressure")
        return Property::entropy_derivative_pressure;
      else if (s == "entropy derivative temperature")
        return Property::entropy_derivative_temperature;
      else if (s == "reaction terms")
        return Property::reaction_terms;
      else
        AssertThrow (false,
                     ExcMessage ("The value <" + s + "> for a material "
                                 "property is not one of the valid values."));
      return Property::viscosity;
    }

    //Copy the requested data for one model
    template <int dim>
    void
    Compositing<dim>::composite(const unsigned int model_index,
                                const typename Interface<dim>::MaterialModelOutputs &evaluated,
                                typename Interface<dim>::MaterialModelOutputs &out) const
    {
      if (model_property_map.at(Property::viscosity) == model_index)
        out.viscosities = evaluated.viscosities;
      if (model_property_map.at(Property::density) == model_index)
        out.densities = evaluated.densities;
      if (model_property_map.at(Property::thermal_expansion_coefficient) == model_index)
        out.thermal_expansion_coefficients = evaluated.thermal_expansion_coefficients;
      if (model_property_map.at(Property::specific_heat) == model_index)
        out.specific_heat = evaluated.specific_heat;
      if (model_property_map.at(Property::compressibility) == model_index)
        out.compressibilities = evaluated.compressibilities;
      if (model_property_map.at(Property::entropy_derivative_pressure) == model_index)
        out.entropy_derivative_pressure = evaluated.entropy_derivative_pressure;
      if (model_property_map.at(Property::entropy_derivative_temperature) == model_index)
        out.entropy_derivative_temperature = evaluated.entropy_derivative_temperature;
      if (model_property_map.at(Property::reaction_terms) == model_index)
        out.reaction_terms == evaluated.reaction_terms;
    }

    template <int dim>
    void
    Compositing<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                               typename Interface<dim>::MaterialModelOutputs &out) const
    {
      typename Interface<dim>::MaterialModelOutputs evaluated(out.viscosities.size(),
                                                              this->introspection().n_compositional_fields);

      for (unsigned int i=0; i<models.size(); ++i)
        {
          models[i]->evaluate(in, evaluated);
          composite(i, evaluated, out);
        }
    }

    template <int dim>
    void
    Compositing<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Compositing");
        {
          prm.declare_entry("Models", "",
                            Patterns::Map(Patterns::Selection("viscosity|density|thermal expansion coefficient|specific heat|compressibility|entropy derivative pressure|entropy derivative temperature|reaction terms"),
                                          Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>())),
                            "The material property and material model associations used for this composite material model");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Compositing<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Compositing");
        {
          model_names.clear();
          const std::vector<std::string> x_models
            = Utilities::split_string_list
              (prm.get("Models"));
          for (std::vector<std::string>::const_iterator p = x_models.begin();
               p != x_models.end(); ++p)
            {
              const std::vector<std::string> split_parts = Utilities::split_string_list(*p, ':');
              AssertThrow (split_parts.size() == 2,
                           ExcMessage("The format for composite models requires that each entry has the form "
                                      "<property> : <model>, but there is no colon in the entry <"
                                      + *p
                                      + ">."));
              Property::MaterialProperty prop = parse_property_name(split_parts[0]);
              AssertThrow(split_parts[1] != "averaging",
                          ExcMessage("You may not use ``averaging'' as the base model for "
                                     "a compositing material model."));
              AssertThrow(split_parts[1] != "compositing",
                          ExcMessage("You may not use ``compositing'' as the base model for "
                                     "a compositing material model."));
              unsigned int model_ind;
              for (model_ind=0; model_ind < model_names.size(); ++model_ind)
                {
                  if (model_names[model_ind] == split_parts[1])
                    break;
                }
              if (model_ind == model_names.size())
                model_names.push_back(split_parts[1]);
              AssertThrow(model_property_map.count(prop)==0,
                          ExcMessage("You may not define a propery multiple times"));

              model_property_map[prop] = model_ind;
            }
          models.resize(model_names.size());

          // create the model and initialize their SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          for (unsigned int i=0; i<model_names.size(); ++i)
            {
              models.at(i).reset(create_material_model<dim>(model_names[i]));
              if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(models.at(i).get()))
                sim->initialize_simulator (this->get_simulator());
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for averaging, it is essential to parse
      parameters related to the base model. */
      for (unsigned int i=0; i<model_names.size(); ++i)
        {
          models[i]->parse_parameters(prm);
          this -> model_dependence.viscosity |= models[i]->get_model_dependence().viscosity;
          this -> model_dependence.density |= models[i]->get_model_dependence().density;
          this -> model_dependence.compressibility |= models[i]->get_model_dependence().compressibility;
          this -> model_dependence.specific_heat |= models[i]->get_model_dependence().specific_heat;
          this -> model_dependence.thermal_conductivity |= models[i]->get_model_dependence().thermal_conductivity;
        }
    }

    template <int dim>
    bool
    Compositing<dim>::
    is_compressible () const
    {
      unsigned int ind = model_property_map.at(Property::compressibility);
      return models[ind]->is_compressible();
    }

    template <int dim>
    double
    Compositing<dim>::
    reference_viscosity() const
    {
      unsigned int ind = model_property_map.at(Property::viscosity);
      return models[ind]->reference_viscosity();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Compositing,
                                   "compositing",
                                   "The ``compositing'' Material model selects material model properties from a "
                                   "given set of other material models."
                                  )
  }
}
