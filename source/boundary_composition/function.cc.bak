/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/boundary_composition/function.h>
#include <aspect/utilities.h>
#include <aspect/global.h>

namespace aspect
{
  namespace BoundaryComposition
  {

    template <int dim>
    double
    Function<dim>::
    boundary_composition (const types::boundary_id /*boundary_indicator*/,
                          const Point<dim> &position,
                          const unsigned int compositional_field) const
    {
      const Utilities::NaturalCoordinate<dim> point =
        this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);
      return function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), compositional_field);
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        function->set_time (this->get_time() / year_in_seconds);
      else
        function->set_time (this->get_time());
    }



    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are 'cartesian', 'spherical', and 'depth'. "
                             "'spherical' coordinates are interpreted as r,phi "
                             "or r,phi,theta in 2d/3d respectively with theta "
                             "being the polar angle. 'depth' will create a "
                             "function, in which only the first parameter is "
                             "non-zero, which is interpreted to be the depth of "
                             "the point.");

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
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = ::aspect::Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));

          try
            {
              function
                = std::make_unique<Functions::ParsedFunction<dim>>(this->n_compositional_fields());
              function->parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Boundary composition model.Function'\n"
                        << "with expression\n"
                        << "\t'" << prm.get("Function expression") << "'"
                        << "More information about the cause of the parse error \n"
                        << "is shown below.\n";
              throw;
            }
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
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(Function,
                                               "function",
                                               "Implementation of a model in which the boundary "
                                               "composition is given in terms of an explicit formula "
                                               "that is elaborated in the parameters in section "
                                               "``Boundary composition model|Function''. "
                                               "\n\n"
                                               "Since the symbol $t$ indicating time "
                                               "may appear in the formulas for the prescribed "
                                               "composition, it is interpreted as having units "
                                               "seconds unless the global input parameter ``Use "
                                               "years in output instead of seconds'' is set, in "
                                               "which case we interpret the formula expressions "
                                               "as having units year."
                                               "\n\n"
                                               "The format of these "
                                               "functions follows the syntax understood by the "
                                               "muparser library, see "
                                               "{ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}
