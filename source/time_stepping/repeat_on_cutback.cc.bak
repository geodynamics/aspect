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
#include <aspect/time_stepping/repeat_on_cutback.h>

namespace aspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    RepeatOnCutback<dim>::execute()
    {
      return std::numeric_limits<double>::max();
    }



    template <int dim>
    std::pair<Reaction, double>
    RepeatOnCutback<dim>::determine_reaction (const TimeStepInfo &info)
    {
      if (info.next_time_step_size < this->get_timestep()*this->repeat_threshold
          && !info.reduced_by_termination_plugin)
        return
          std::make_pair<Reaction, double>(Reaction::repeat_step,
                                           this->get_timestep()*this->cut_back_amount);
      else
        return
          std::make_pair<Reaction, double>(Reaction::advance,
                                           std::numeric_limits<double>::max());
    }



    template <int dim>
    void
    RepeatOnCutback<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Repeat on cutback");

        prm.declare_entry("Relative repeat threshold", "0.2",
                          Patterns::Double (0.),
                          "A factor that controls when a step is going to be repeated. If "
                          "the newly computed step size is smaller than the last step size "
                          "multiplied by this factor, the step is repeated.");

        prm.declare_entry("Cut back amount", "0.5",
                          Patterns::Double (0.),
                          "A factor that controls the size of the time step when repeating. The "
                          "default of 0.5 corresponds to 50\\% of the original step taken.");

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    RepeatOnCutback<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Repeat on cutback");

        repeat_threshold = prm.get_double("Relative repeat threshold");
        cut_back_amount = prm.get_double("Cut back amount");

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
    ASPECT_REGISTER_TIME_STEPPING_MODEL(RepeatOnCutback,
                                        "repeat on cutback",
                                        "This time stepping plugin will detect a situation where the computed time step "
                                        "shrinks by more than a user-controlled factor. In that situation, the previous "
                                        "time step will be repeated with a smaller step size.\n"
                                        "A large reduction in time step size typically happens when velocities change "
                                        "abruptly. Repeating the time step ensure properly resolving this event. It is "
                                        "useful to consider setting the \"Maximum relative increase in time step\" option "
                                        "to avoid repeatedly repeating every other time step.")
  }
}
