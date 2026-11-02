// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/boundary_temperature/initial_temperature.h>
#include <aspect/initial_temperature/interface.h>


namespace aspect
{
  namespace BoundaryTemperature
  {
// ------------------------------ InitialTemperature -------------------

    template <int dim>
    void
    InitialTemperature<dim>::initialize()
    {
      // Make sure we keep track of the initial temperature manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_temperature = this->get_initial_temperature_manager_pointer();
    }



    template <int dim>
    double
    InitialTemperature<dim>::
    boundary_temperature (const types::boundary_id,
                          const Point<dim> &position) const
    {
      return initial_temperature->initial_temperature(position);
    }


    template <int dim>
    double
    InitialTemperature<dim>::
    minimal_temperature (const std::set<types::boundary_id> &) const
    {
      return min_temperature;
    }



    template <int dim>
    double
    InitialTemperature<dim>::
    maximal_temperature (const std::set<types::boundary_id> &) const
    {
      return max_temperature;
    }



    template <int dim>
    void
    InitialTemperature<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Initial temperature");
        {
          prm.declare_entry ("Minimal temperature", "0.",
                             Patterns::Double (),
                             "Minimal temperature. Units: $\\text{K}$.");
          prm.declare_entry ("Maximal temperature", "3773.",
                             Patterns::Double (),
                             "Maximal temperature. Units: $\\text{K}$.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    InitialTemperature<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Initial temperature");
        {
          min_temperature = prm.get_double ("Minimal temperature");
          max_temperature = prm.get_double ("Maximal temperature");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(InitialTemperature,
                                               "initial temperature",
                                               "A model in which the temperature at the boundary "
                                               "is chosen to be the same as given in the initial "
                                               "conditions."
                                               "\n\n"
                                               "Because this class simply takes what the initial "
                                               "temperature had described, this class can not "
                                               "know certain pieces of information such as the "
                                               "minimal and maximal temperature on the boundary. "
                                               "For operations that require this, for example in "
                                               "post-processing, this boundary temperature model "
                                               "must therefore be told what the minimal and "
                                               "maximal values on the boundary are. This is done "
                                               "using parameters set in section ``Boundary temperature model/Initial temperature''.")
  }
}
