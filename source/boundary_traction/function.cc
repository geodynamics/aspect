/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/boundary_traction/function.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace BoundaryTraction
  {
    template <int dim>
    Function<dim>::Function ()
      :
      boundary_traction_function (dim)
    {}



    template <int dim>
    Tensor<1,dim>
    Function<dim>::
    boundary_traction (const types::boundary_id,
                       const Point<dim> &p,
                       const Tensor<1,dim> &) const
    {
      Tensor<1,dim> traction;
      if (coordinate_system == Utilities::Coordinates::cartesian)
        {
          for (unsigned int d=0; d<dim; ++d)
            traction[d] = boundary_traction_function.value(p,d);

        }
      else if (coordinate_system == Utilities::Coordinates::spherical)
        {
          const std_cxx11::array<double,dim> spherical_coordinates =
            aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(p);
          Point<dim> point;

          for (unsigned int d=0; d<dim; ++d)
            point[d] = spherical_coordinates[d];

          for (unsigned int d=0; d<dim; ++d)
            traction[d] = boundary_traction_function.value(point,d);

        }
      else if (coordinate_system == Utilities::Coordinates::depth)
        {
          const double depth = this->get_geometry_model().depth(p);
          Point<dim> point;
          point(0) = depth;

          for (unsigned int d=0; d<dim; ++d)
            traction[d] = boundary_traction_function.value(point,d);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
          return numbers::signaling_nan<Tensor<1,dim> >();
        }

      return traction;
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_traction_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_traction_function.set_time (this->get_time());
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }
        try
          {
            boundary_traction_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary traction model.Function'\n"
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
  namespace BoundaryTraction
  {
    ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(Function,
                                            "function",
                                            "Implementation of a model in which the boundary "
                                            "traction is given in terms of an explicit formula "
                                            "that is elaborated in the parameters in section "
                                            "``Boundary traction model|Function''. "
                                            "\n\n"
                                            "The formula you describe in the mentioned "
                                            "section is a semicolon separated list of traction components "
                                            "for each of the $d$ components of the traction vector. "
                                            "These $d$ formulas are interpreted as having units "
                                            "Pa.")
  }
}
