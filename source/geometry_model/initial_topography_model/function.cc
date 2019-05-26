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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/geometry_model/initial_topography_model/function.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace InitialTopographyModel
  {
    template <int dim>
    Function<dim>::Function ()
      :
    max_topo(0),
      initial_topography_function (1),
      coordinate_system(Utilities::Coordinates::CoordinateSystem::cartesian)
    {}

    template <int dim>
    double
    Function<dim>::
    value (const Point<dim-1> &surface_point) const
    {
      Point<dim> global_point;
      if (dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()) != nullptr)
        {
          // No need to set the vertical coordinate correctly,
          // because it will be thrown away in get_data_component anyway
          for (unsigned int d=0; d<dim-1; d++)
            global_point[d] = surface_point[d];
        }
      else if (dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != nullptr ||
               dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != nullptr ||
               dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != nullptr)
        {
          // No need to set the radial coordinate correctly,
          // because it will be thrown away in get_data_component anyway
          std::array<double, dim> point;
          point[0] = 6371000.0;
          for (unsigned int d=0; d<dim-1; d++)
            point[d+1] = surface_point[d];

          global_point = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(point);
        }
      else
        AssertThrow(false, ExcNotImplemented());

//      Utilities::NaturalCoordinate<dim> point =
//        this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);
//      const double topo = initial_topography_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()));
      const double topo = initial_topography_function.value(global_point);

      return topo;
    }


    template <int dim>
    double
    Function<dim>::
    max_topography () const
    {
      return 0;
    }

    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
      prm.enter_subsection("Initial topography model");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry ("Maximum topography value", "2000",
                             Patterns::Double (0),
                             "The maximum value the topography given by "
                             "the function can take. ");
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle. "); 

          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
      prm.enter_subsection("Initial topography model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }
        max_topo = prm.get_double("Maximum topography value");
        try
          {
            initial_topography_function.parse_parameters (prm);
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
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(Function,
                                             "function",
                                             "Implementation of a model in which the initial topography "
                                             "is described by a function. ")
  }
}
