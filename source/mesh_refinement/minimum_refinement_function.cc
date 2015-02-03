/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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



#include <aspect/mesh_refinement/minimum_refinement_function.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <math.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    MinimumRefinementFunction<dim>::tag_additional_cells () const
    {
      for (typename Triangulation<dim>::active_cell_iterator
           cell = this->get_triangulation().begin_active();
           cell != this->get_triangulation().end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
              bool refine = false;
              bool clear_coarsen = false;

              for ( unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;  ++v)
                {
                  const Point<dim> vertex = cell->vertex(v);
                  double minimum_refinement_level = 0;

                  if (coordinate_system == depth)
                    {
                      const double depth = this->get_geometry_model().depth(vertex);
                      Point<dim> point;
                      point(0) = depth;
                      minimum_refinement_level = min_refinement_level.value(point);
                    }
                  else if (coordinate_system == spherical)
                    {
                      const std_cxx1x::array<double,dim> spherical_coordinates =
                        aspect::Utilities::spherical_coordinates(vertex);

                      // Conversion to evaluate the spherical coordinates in the minimum
                      // refinement level function.
                      Point<dim> point;
                      for (unsigned int i = 0; i<dim; ++i)
                        point[i] = spherical_coordinates[i];

                      minimum_refinement_level = min_refinement_level.value(point);
                    }
                  else if (coordinate_system == cartesian)
                    {
                      minimum_refinement_level = min_refinement_level.value(vertex);
                    }

                  if (cell->level() <= rint(minimum_refinement_level))
                    clear_coarsen = true;
                  if (cell->level() <  rint(minimum_refinement_level))
                    {
                      refine = true;
                      break;
                    }
                }

              if (clear_coarsen)
                cell->clear_coarsen_flag ();
              if (refine)
                cell->set_refine_flag ();
            }
        }
    }

    template <int dim>
    void
    MinimumRefinementFunction<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {

        prm.enter_subsection("Minimum refinement function");
        {
          /**
           * Choose the coordinates to evaluate the minimum refinement level
           * function. The function can be declared in dependence of depth,
           * cartesian coordinates or spherical coordinates. Note that the order
           * of spherical coordinates is r,phi,theta and not r,theta,phi, since
           * this allows for dimension independent expressions.
           */
          prm.declare_entry ("Coordinate system", "depth",
                             Patterns::Selection ("depth|cartesian|spherical"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are 'depth', 'cartesian' and 'spherical'. 'depth' "
                             "will create a function, in which only the first "
                             "variable is non-zero, which is interpreted to "
                             "be the depth of the point. 'spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle.");
          /**
           * Let the function that describes the minimal level of refinement
           * as a function of position declare its parameters.
           * This defines the minimum refinement level each cell should have,
           * and that can not be exceeded by coarsening.
           */
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MinimumRefinementFunction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Minimum refinement function");
        {
          if (prm.get ("Coordinate system") == "depth")
            coordinate_system = depth;
          else if (prm.get ("Coordinate system") == "cartesian")
            coordinate_system = cartesian;
          else if (prm.get ("Coordinate system") == "spherical")
            coordinate_system = spherical;
          else
            AssertThrow (false, ExcNotImplemented());

          try
            {
              min_refinement_level.parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Mesh refinement.Minimum refinement function'\n"
                        << "with expression\n"
                        << "\t'" << prm.get("Function expression") << "'";
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
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(MinimumRefinementFunction,
                                              "minimum refinement function",
                                              "A mesh refinement criterion that ensures a "
                                              "minimum refinement level described by an "
                                              "explicit formula with the depth or position "
                                              "as argument. Which coordinate representation "
                                              "is used is determined by an input parameter. "
                                              "Whatever the coordinate system chosen, the "
                                              "function you provide in the input file will "
                                              "by default depend on variables 'x', 'y' and "
                                              "'z' (if in 3d). However, the meaning of these "
                                              "symbols depends on the coordinate system. In "
                                              "the Cartesian coordinate system, they simply "
                                              "refer to their natural meaning. If you have "
                                              "selected 'depth' for the coordinate system, "
                                              "then 'x' refers to the depth variable and 'y' "
                                              "and 'z' will simply always be zero. If you "
                                              "have selected a spherical coordinate system, "
                                              "then 'x' will refer to the radial distance of "
                                              "the point to the origin, 'y' to the azimuth "
                                              "angle and 'z' to the polar angle measured "
                                              "positive from the north pole. Note that the "
                                              "order of spherical coordinates is r,phi,theta "
                                              "and not r,theta,phi, since this allows for "
                                              "dimension independent expressions. "
                                              "After evaluating the function, its values are "
                                              "rounded to the nearest integer."
                                              "\n\n"
                                              "The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
