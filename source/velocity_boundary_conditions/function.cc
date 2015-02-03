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


#include <aspect/velocity_boundary_conditions/function.h>
#include <aspect/global.h>

namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    template <int dim>
    Function<dim>::Function ()
      :
      boundary_velocity_function (dim)
    {}



    template <int dim>
    Tensor<1,dim>
    Function<dim>::
    boundary_velocity (const Point<dim> &p) const
    {
      Tensor<1,dim> velocity;
      for (unsigned int d=0; d<dim; ++d)
        velocity[d] = boundary_velocity_function.value(p,d);

      // Aspect always wants things in MKS system. however, as described
      // in the documentation of this class, we interpret the formulas
      // given to this plugin as meters per year if the global flag
      // for using years instead of seconds is given. so if someone
      // write "5" in their parameter file and sets the flag, then this
      // means "5 meters/year" and we need to convert it to the Aspect
      // time system by dividing by the number of seconds per year
      // to get MKS units
      if (this->convert_output_to_years())
        return velocity / year_in_seconds;
      else
        return velocity;
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_velocity_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_velocity_function.set_time (this->get_time());
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("Function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("Function");
        try
          {
            boundary_velocity_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary velocity model.Function'\n"
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
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(Function,
                                                 "function",
                                                 "Implementation of a model in which the boundary "
                                                 "velocity is given in terms of an explicit formula "
                                                 "that is elaborated in the parameters in section "
                                                 "``Boundary velocity model|Function''. The format of these "
                                                 "functions follows the syntax understood by the "
                                                 "muparser library, see Section~\\ref{sec:muparser-format}."
                                                 "\n\n"
                                                 "The formula you describe in the mentioned "
                                                 "section is a semicolon separated list of velocities "
                                                 "for each of the $d$ component of the velocity vector. "
                                                 "These $d$ formulas are interpreted as having units "
                                                 "m/s, unless the global input parameter ``Use "
                                                 "years in output instead of seconds'' is set, in "
                                                 "which case we interpret the formula expressions "
                                                 "as having units m/year."
                                                 "\n\n"
                                                 "Likewise, since the symbol $t$ indicating time "
                                                 "may appear in the formulas for the prescribed "
                                                 "velocities, it is interpreted as having units "
                                                 "seconds unless the global parameter above has "
                                                 "been set.")
  }
}
