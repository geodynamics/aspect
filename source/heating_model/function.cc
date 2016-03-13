/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/heating_model/function.h>
#include <aspect/global.h>

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
          // return a constant value
          const Point<dim> position = material_model_inputs.position[q];
          heating_model_outputs.heating_source_terms[q] = heating_model_function.value(position)
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
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Function");
        {
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
        {
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
                                  "muparser library, see Section~\\ref{sec:muparser-format}."
                                  "\n\n"
                                  "The formula is interpreted as having units "
                                  "W/kg."
                                  "\n\n"
                                  "Since the symbol $t$ indicating time "
                                  "may appear in the formulas for the heating rate"
                                  ", it is interpreted as having units "
                                  "seconds unless the global parameter ``Use "
                                  "years in output instead of seconds'' is set.")
  }
}
