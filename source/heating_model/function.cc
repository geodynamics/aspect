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


#include <aspect/heating_model/function.h>

#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    Function<dim>::Function ()
      :
      heating_model_function (1)
    {}



    template <int dim>
    void
    Function<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // convert the position into the selected coordinate system
          const Point<dim> position = material_model_inputs.position[q];
          const Utilities::NaturalCoordinate<dim> point =
            this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);

          // then compute the heating function value at this position
          heating_model_outputs.heating_source_terms[q] = heating_model_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()))
                                                          * material_model_outputs.densities[q];
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    Function<dim>::update ()
    {
      const double time = this->get_time();
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        heating_model_function.set_time (time / year_in_seconds);
      else
        heating_model_function.set_time (time);
    }



    template <int dim>
    MaterialModel::MaterialProperties::Property
    Function<dim>::get_required_properties () const
    {
      return MaterialModel::MaterialProperties::none;
    }



    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2d/3d "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Function");
        {
          try
            {
              heating_model_function.parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Heating model.Function'\n"
                        << "with expression\n"
                        << "\t'" << prm.get("Function expression") << "'"
                        << "More information about the cause of the parse error \n"
                        << "is shown below.\n";
              throw;
            }

          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
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
    ASPECT_REGISTER_HEATING_MODEL(Function,
                                  "function",
                                  "Implementation of a model in which the heating "
                                  "rate is given in terms of an explicit formula "
                                  "that is elaborated in the parameters in section "
                                  "``Heating model|Function''. The format of these "
                                  "functions follows the syntax understood by the "
                                  "muparser library, see "
                                  "{ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`."
                                  "\n\n"
                                  "The formula is interpreted as having units "
                                  "W/kg."
                                  "\n\n"
                                  "Since the symbol $t$ indicating time "
                                  "may appear in the formulas for the heating rate"
                                  ", it is interpreted as having units "
                                  "seconds unless the global parameter "
                                  "``Use years instead of seconds'' is set.")
  }
}
