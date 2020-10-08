/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

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
      double new_time_step = function.value(Point<1>(this->get_time()));

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
                          "Expression for the time step size as a function of 'time'.");
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
                      << "\t'Heating model.Function'\n"
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
