// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/termination_criteria/end_walltime.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    // The start_walltime is made as static and initialized here to
    // make sure it is initialized when the program starts.
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
                           "The wall time of the simulation. Unit: \\si{\\hour}.");
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
