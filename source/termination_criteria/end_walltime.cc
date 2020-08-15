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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/termination_criteria/end_walltime.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    // The start_walltime is made as static and initialized here to
    // make sure it is initialized right way as the program starts.
    template <int dim>
    std::time_t
    EndWalltime<dim>::start_walltime = std::time(nullptr);

    template <int dim>
    bool
    EndWalltime<dim>::execute()
    {
      return std::difftime(std::time(nullptr),start_walltime) >= walltime_duration;
    }


    template <int dim>
    void
    EndWalltime<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.declare_entry ("Wall time",
                           "24.",
                           Patterns::Double (0.),
                           "The wall time of the simulation. Unit: hours.");
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    EndWalltime<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        // Change from hours to seconds:
        walltime_duration = static_cast<unsigned int>(prm.get_double ("Wall time") * 3600.);
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
    ASPECT_REGISTER_TERMINATION_CRITERION(EndWalltime,
                                          "wall time",
                                          "Terminate the simulation once the wall time limit has reached.")
  }
}
