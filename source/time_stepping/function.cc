// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/time_stepping/function.h>

namespace aspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    Function<dim>::execute()
    {
      const double time_scale    = this->convert_output_to_years() ? constants::year_in_seconds : 1.0;
      const double new_time_step = function.value(Point<1>(this->get_time()) / time_scale) * time_scale;

      AssertThrow (new_time_step > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length of the \"function\" plugin was: " +
                              std::to_string(new_time_step) + "."));

      return new_time_step;
    }



    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Function");
        Functions::ParsedFunction<1>::declare_parameters (prm, 1);
        prm.declare_entry("Function expression", "1.0",
                          Patterns::Anything(),
                          "Expression for the time step size as a function of 'time'."
                          "Units: \\si{\\year} if the "
                          "'Use years instead of seconds' parameter is set; "
                          "\\si{\\second} otherwise.");
        prm.declare_entry("Variable names", "time",
                          Patterns::Anything(),
                          "Name for the variable representing the current time.");
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Function");
        try
          {
            function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Timestepping.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
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
  namespace TimeStepping
  {
    ASPECT_REGISTER_TIME_STEPPING_MODEL(Function,
                                        "function",
                                        "This model uses a time step specified in the parameter "
                                        "file specified as a function of time. This plugin will always "
                                        "request advancing to the next time step.")
  }
}
