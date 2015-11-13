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


#include <aspect/global.h>
#include <aspect/prescribed_stokes_solution/function.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    template <int dim>
    Function<dim>::Function ()
      :
      prescribed_velocity_function (dim),
      prescribed_pressure_function (1)
    {}

    template <int dim>
    void
    Function<dim>::stokes_solution (const Point<dim> &position, Vector<double> &value) const
    {
      //velocity
      for (unsigned int d=0; d<dim; ++d)
        value[d] = prescribed_velocity_function.value(position,d);

      // pressure
      value(dim) = prescribed_pressure_function.value(position);
    }

    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        {
          prescribed_velocity_function.set_time (this->get_time() / year_in_seconds);
          prescribed_pressure_function.set_time (this->get_time() / year_in_seconds);
        }
      else
        {
          prescribed_velocity_function.set_time (this->get_time());
          prescribed_pressure_function.set_time (this->get_time());
        }
    }

    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed Stokes solution");
      {
        prm.enter_subsection("Velocity function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
        prm.enter_subsection("Pressure function");
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
      prm.enter_subsection("Prescribed Stokes solution");
      {
        prm.enter_subsection("Velocity function");
        try
          {
            prescribed_velocity_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Initial conditions.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
          }
        prm.leave_subsection();
        prm.enter_subsection("Pressure function");
        try
          {
            prescribed_pressure_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Initial conditions.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
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
  namespace PrescribedStokesSolution
  {
    ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION(Function,
                                               "function",
                                               "This plugin allows to prescribe the Stokes solution "
                                               "for the velocity and pressure field in terms of an "
                                               "explicit formula. The format of these "
                                               "functions follows the syntax understood by the "
                                               "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}

