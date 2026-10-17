// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2017 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/adiabatic_conditions/function.h>

namespace aspect
{
  namespace AdiabaticConditions
  {
    template <int dim>
    Function<dim>::Function()
      : function (3)
    {}

    template <int dim>
    void
    Function<dim>::initialize()
    {
    }

    template <int dim>
    bool
    Function<dim>::is_initialized() const
    {
      return true;
    }



    template <int dim>
    double Function<dim>::pressure (const Point<dim> &p) const
    {
      const double z = this->get_geometry_model().depth(p);
      return function.value(Point<1>(z), 1);
    }



    template <int dim>
    double Function<dim>::temperature (const Point<dim> &p) const
    {
      const double z = this->get_geometry_model().depth(p);
      return function.value(Point<1>(z), 0);
    }

    template <int dim>
    double Function<dim>::density (const Point<dim> &p) const
    {
      const double z = this->get_geometry_model().depth(p);
      return function.value(Point<1>(z), 2);
    }



    template <int dim>
    double Function<dim>::density_derivative (const Point<dim> &p) const
    {
      // TODO: better eps or make it a user input
      const double z = this->get_geometry_model().depth(p);
      const double z2 = z + (1.e6 * std::numeric_limits<double>::epsilon())
                        * this->get_geometry_model().maximal_depth();
      return (function.value(Point<1>(z), 2)
              - function.value(Point<1>(z2), 2))/(z-z2);
    }

    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Function");
        Functions::ParsedFunction<1>::declare_parameters (prm, 3);
        prm.declare_entry("Function expression","0.0; 0.0; 1.0",
                          Patterns::Anything(),
                          "Expression for the adiabatic temperature, "
                          "pressure, and density separated by "
                          "semicolons as a function of `depth'.");
        prm.declare_entry("Variable names","depth");
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Function");
        try
          {
            function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Adiabatic conditions model.Function'\n"
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
  namespace AdiabaticConditions
  {
    ASPECT_REGISTER_ADIABATIC_CONDITIONS_MODEL(Function,
                                               "function",
                                               "A model in which the adiabatic profile is "
                                               "specified by a user defined function. The "
                                               "supplied function has to contain "
                                               "temperature, pressure, and density as a function "
                                               "of depth in this order.")
  }
}
