/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/geometry_model/two_merged_boxes.h>
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
      // In a first step, create a global 'dim'-dimensional point that we can pass to the
      // function expression as input -- because the function is a dim-dimensional
      // function.
      //
      // Different geometry models pass the surface point in in different ways.
      // In the following, we will first normalize the input to a dim-dimensional
      // point with a dummy vertical/radial coordinate that, one would hope,
      // the function expression will then simply ignore.
      Point<dim> cartesian_point;

      if (Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()) ||
          Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
        {
          for (unsigned int d=0; d<dim-1; ++d)
            cartesian_point[d] = surface_point[d];

          // Now for the vertical component:
          cartesian_point[dim-1] = 0;

          // The point as it is would have to be translated into a different
          // coordinate system if that was requested in the input file.
          Assert (coordinate_system == Utilities::Coordinates::CoordinateSystem::cartesian,
                  ExcNotImplemented());
        }
      else if (Plugins::plugin_type_matches<GeometryModel::Sphere<dim>>(this->get_geometry_model()) ||
               Plugins::plugin_type_matches<GeometryModel::SphericalShell<dim>>(this->get_geometry_model()) ||
               Plugins::plugin_type_matches<GeometryModel::Chunk<dim>>(this->get_geometry_model()) )
        {
          std::array<double, dim> spherical_point;

          // This value does not affect the actual topography.
          spherical_point[0] = 6371000.0;
          for (unsigned int d=0; d<dim-1; ++d)
            spherical_point[d+1] = surface_point[d];

          cartesian_point = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(spherical_point);
        }
      else
        AssertThrow(false, ExcNotImplemented());

      Point<dim> function_point;
      if (coordinate_system == Utilities::Coordinates::CoordinateSystem::cartesian)
        {
          function_point = cartesian_point;
        }
      else if (coordinate_system == Utilities::Coordinates::CoordinateSystem::spherical)
        {
          const std::array<double, dim> spherical_coords =
            Utilities::Coordinates::cartesian_to_spherical_coordinates(cartesian_point);
          function_point = Utilities::convert_array_to_point<dim>(spherical_coords);
        }
      else
        AssertThrow(false, ExcMessage("Unsupported coordinate system in initial topography function"));

      const double topo = initial_topography_function.value(function_point);

      return topo;
    }


    template <int dim>
    double
    Function<dim>::
    max_topography () const
    {
      return max_topo;
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
            prm.declare_entry ("Maximum topography value", "2000.",
                               Patterns::Double (0.),
                               "The maximum value the topography given by "
                               "the function can take. ");
            prm.declare_entry ("Coordinate system", "cartesian",
                               Patterns::Selection ("cartesian|spherical"),
                               "A selection that determines the assumed coordinate "
                               "system for the function variables. Allowed values "
                               "are `cartesian' and `spherical'. `spherical' coordinates "
                               "are interpreted as r,phi or r,phi,theta in 2d/3d "
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
                                             "is described by a function in Cartesian or spherical coordinates. ")
  }
}
