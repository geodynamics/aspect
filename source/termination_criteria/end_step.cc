// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2018 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
