/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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



#include <aspect/mesh_refinement/maximum_refinement_function.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <cmath>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    MaximumRefinementFunction<dim>::update ()
    {
      const double time = this->get_time() /
                          (this->convert_output_to_years()
                           ?
                           year_in_seconds
                           :
                           1.0);

      max_refinement_level.set_time(time);
    }

    template <int dim>
    void
    MaximumRefinementFunction<dim>::tag_additional_cells () const
    {
      for (const auto &cell : this->get_triangulation().active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              bool coarsen = false;
              bool clear_refine = false;

              for ( unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;  ++v)
                {
                  const Point<dim> vertex = cell->vertex(v);
                  Utilities::NaturalCoordinate<dim> point =
                    this->get_geometry_model().cartesian_to_other_coordinates(vertex, coordinate_system);

                  const double maximum_refinement_level = max_refinement_level.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()));

                  if (cell->level() >= rint(maximum_refinement_level))
                    clear_refine = true;
                  if (cell->level() >  rint(maximum_refinement_level))
                    {
                      coarsen = true;
                      break;
                    }
                }

              if (clear_refine)
                cell->clear_refine_flag ();
              if (coarsen)
                cell->set_coarsen_flag ();
            }
        }
    }

    template <int dim>
    void
    MaximumRefinementFunction<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {

        prm.enter_subsection("Maximum refinement function");
        {
          /**
           * Choose the coordinates to evaluate the maximum refinement level
           * function. The function can be declared in dependence of depth,
           * cartesian coordinates or spherical coordinates. Note that the order
           * of spherical coordinates is r,phi,theta and not r,theta,phi, since
           * this allows for dimension independent expressions.
           */
          prm.declare_entry ("Coordinate system", "depth",
                             Patterns::Selection ("depth|cartesian|spherical"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `depth', `cartesian' and `spherical'. `depth' "
                             "will create a function, in which only the first "
                             "variable is non-zero, which is interpreted to "
                             "be the depth of the point. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle.");
          /**
           * Let the function that describes the maximal level of refinement
           * as a function of position declare its parameters.
           * This defines the maximum refinement level each cell should have,
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
    MaximumRefinementFunction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Maximum refinement function");
        {
          if (prm.get ("Coordinate system") == "depth")
            coordinate_system = Utilities::Coordinates::depth;
          else if (prm.get ("Coordinate system") == "cartesian")
            coordinate_system = Utilities::Coordinates::cartesian;
          else if (prm.get ("Coordinate system") == "spherical")
            coordinate_system = Utilities::Coordinates::spherical;
          else
            AssertThrow (false, ExcNotImplemented());

          try
            {
              max_refinement_level.parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Mesh refinement.Maximum refinement function'\n"
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
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(MaximumRefinementFunction,
                                              "maximum refinement function",
                                              "A mesh refinement criterion that ensures a "
                                              "maximum refinement level described by an "
                                              "explicit formula with the depth or position "
                                              "as argument. Which coordinate representation "
                                              "is used is determined by an input parameter. "
                                              "Whatever the coordinate system chosen, the "
                                              "function you provide in the input file will "
                                              "by default depend on variables `x', `y' and "
                                              "`z' (if in 3d). However, the meaning of these "
                                              "symbols depends on the coordinate system. In "
                                              "the Cartesian coordinate system, they simply "
                                              "refer to their natural meaning. If you have "
                                              "selected `depth' for the coordinate system, "
                                              "then `x' refers to the depth variable and `y' "
                                              "and `z' will simply always be zero. If you "
                                              "have selected a spherical coordinate system, "
                                              "then `x' will refer to the radial distance of "
                                              "the point to the origin, `y' to the azimuth "
                                              "angle and `z' to the polar angle measured "
                                              "positive from the north pole. Note that the "
                                              "order of spherical coordinates is r,phi,theta "
                                              "and not r,theta,phi, since this allows for "
                                              "dimension independent expressions. "
                                              "Each coordinate system also includes a final `t' "
                                              "variable which represents the model time, evaluated "
                                              "in years if the 'Use years in output instead of seconds' "
                                              "parameter is set, otherwise evaluated in seconds. "
                                              "After evaluating the function, its values are "
                                              "rounded to the nearest integer."
                                              "\n\n"
                                              "The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
