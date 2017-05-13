/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/function.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <deal.II/base/signaling_nan.h>

#include <iostream>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    Function<dim>::Function ()
      :
      function (1)
    {}

    template <int dim>
    double
    Function<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      if (coordinate_system == ::aspect::Utilities::Coordinates::CoordinateSystem::cartesian)
        {
          return function.value(position);
        }
      else if (coordinate_system == ::aspect::Utilities::Coordinates::CoordinateSystem::spherical)
        {
          const std_cxx11::array<double,dim> spherical_coordinates =
            aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
          Point<dim> point;

          for (unsigned int i = 0; i<dim; ++i)
            point[i] = spherical_coordinates[i];

          return function.value(point);
        }
          else if (coordinate_system == ::aspect::Utilities::Coordinates::CoordinateSystem::depth)
        {
          const double depth = this->get_geometry_model().depth(position);
          Point<dim> point;
          point(0) = depth;
          return function.value(point);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
          return numbers::signaling_nan<double>();
        }
    }

    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Function");
        {
          /**
           * Choose the coordinates to evaluate the maximum refinement level
           * function. The function can be declared in dependence of depth,
           * cartesian coordinates or spherical coordinates. Note that the order
           * of spherical coordinates is r,phi,theta and not r,theta,phi, since
           * this allows for dimension independent expressions.
           */
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are 'cartesian', 'spherical', and 'depth'. 'spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle.'depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

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
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = ::aspect::Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }

        try
          {
            function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Initial conditions.Function'\n"
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
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Function,
                                              "function",
                                              "Specify the initial temperature in terms of an "
                                              "explicit formula. The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
