// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2023-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/time_stepping/repeat_on_nonlinear_fail.h>

namespace aspect
{
  namespace TimeStepping
  {
    template <int dim>
    RepeatOnNonlinearFail<dim>::RepeatOnNonlinearFail()
      : nonlinear_solver_just_failed (false),
        maximum_number_of_repeats (10),
        current_number_of_repeats (0)
    {}



    template <int dim>
    double
    RepeatOnNonlinearFail<dim>::execute()
    {
      return std::numeric_limits<double>::max();
    }



    template <int dim>
    void
    RepeatOnNonlinearFail<dim>::nonlinear_solver_has_failed() const
    {
      nonlinear_solver_just_failed = true;
    }



    template <int dim>
    std::pair<Reaction, double>
    RepeatOnNonlinearFail<dim>::determine_reaction (const TimeStepInfo &/*info*/)
    {
      if (nonlinear_solver_just_failed)
        {
          nonlinear_solver_just_failed = false;
          ++current_number_of_repeats;

          AssertThrow(current_number_of_repeats < maximum_number_of_repeats,
                      ExcMessage("The nonlinear solver did not converge even after cutting the timestep "
                                 "size several times."));

          return
            std::make_pair<Reaction, double>(Reaction::repeat_step,
                                             this->get_timestep()*this->cut_back_factor);
        }
      else
        {
          current_number_of_repeats = 0;
          return
            std::make_pair<Reaction, double>(Reaction::advance,
                                             std::numeric_limits<double>::max());
        }
    }



    template <int dim>
    void
    RepeatOnNonlinearFail<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Repeat on nonlinear solver failure");

        prm.declare_entry("Cut back factor", "0.5",
                          Patterns::Double (0.),
                          "A factor that controls the size of the time step when repeating. The "
                          "default of 0.5 corresponds to 50\\% of the original step taken.");

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    RepeatOnNonlinearFail<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Repeat on nonlinear solver failure");

        cut_back_factor = prm.get_double("Cut back factor");

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
    ASPECT_REGISTER_TIME_STEPPING_MODEL(RepeatOnNonlinearFail,
                                        "repeat on nonlinear solver failure",
                                        "This time stepping plugin will react when the nonlinear solver "
                                        "does not converge in the specified maximum number of iterations "
                                        "and repeats the current timestep with a smaller step size. "
                                        "This plugin is enabled automatically if \"Nonlinear solver "
                                        "failure strategy\" is set to \"cut timestep size\".")
  }
}
