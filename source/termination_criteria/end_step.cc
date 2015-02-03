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

#include <aspect/termination_criteria/end_step.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    EndStep<dim>::execute()
    {
      return (this->get_timestep_number () > end_step);
    }

    template <int dim>
    void
    EndStep<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.declare_entry ("End step", "100",
                           Patterns::Integer (0),
                           "Terminate the simulation once the specified timestep has been reached.");
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    EndStep<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        end_step = prm.get_integer ("End step");
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TerminationCriteria
  {
    ASPECT_REGISTER_TERMINATION_CRITERION(EndStep,
                                          "end step",
                                          "Terminate the simulation once the specified timestep "
                                          "has been reached. ")
  }
}
