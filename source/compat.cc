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


#include <aspect/global.h>

#include <aspect/compat.h>

// deal.II versions up to 9.6 had a bug for very thin shell geometries.
// This function contains a fixed version.
#if !DEAL_II_VERSION_GTE(9, 7, 0)

#  include <deal.II/grid/grid_generator.h>

namespace aspect
{
  template <int dim>
  void
  colorize_quarter_hyper_shell(Triangulation<dim> & /*tria*/,
                               const Point<dim> & /*center*/,
                               const double /*inner_radius*/,
                               const double /*outer_radius*/)
  {
    AssertThrow(false, ExcNotImplemented());
  }

  /**
   * Assign boundary number zero the inner shell boundary, one to the outer
   * shell boundary, two to the face with x=0, three to the face with y=0,
   * four to the face with z=0.
   */
  template <>
  void
  colorize_quarter_hyper_shell(Triangulation<3> &tria, const Point<3> &center, const double inner_radius, const double outer_radius)
  {
    if (tria.n_cells() != 3)
      AssertThrow(false, ExcNotImplemented());

    const double middle_radius = (outer_radius - inner_radius) / 2e0 + inner_radius;
    const double eps           = 1e-3 * middle_radius;

    for (const auto &cell : tria.cell_iterators())
      for (const unsigned int f : cell->face_indices())
        {
          const auto face = cell->face(f);
          if (!face->at_boundary())
            continue;

          const double face_center_radius = (face->center(true) - center).norm();

          if (std::fabs(face->center()[0]) < eps) // x = 0 set boundary 2
            {
              face->set_boundary_id(2);
              for (const unsigned int line_no : face->line_indices())
                if (face->line(line_no)->at_boundary())
                  if (std::fabs(face->line(line_no)->vertex(0).norm() - face->line(line_no)->vertex(1).norm()) > eps)
                    face->line(line_no)->set_boundary_id(2);
            }
          else if (std::fabs(face->center()[1]) < eps) // y = 0 set boundary 3
            {
              face->set_boundary_id(3);
              for (const unsigned int line_no : face->line_indices())
                if (face->line(line_no)->at_boundary())
                  if (std::fabs(face->line(line_no)->vertex(0).norm() - face->line(line_no)->vertex(1).norm()) > eps)
                    face->line(line_no)->set_boundary_id(3);
            }
          else if (std::fabs(face->center()[2]) < eps) // z = 0 set boundary 4
            {
              face->set_boundary_id(4);
              for (const unsigned int line_no : face->line_indices())
                if (face->line(line_no)->at_boundary())
                  if (std::fabs(face->line(line_no)->vertex(0).norm() - face->line(line_no)->vertex(1).norm()) > eps)
                    face->line(line_no)->set_boundary_id(4);
            }
          else if (face_center_radius < middle_radius) // inner radius set boundary 0
            {
              face->set_boundary_id(0);
              for (const unsigned int line_no : face->line_indices())
                if (face->line(line_no)->at_boundary())
                  if (std::fabs(face->line(line_no)->vertex(0).norm() - face->line(line_no)->vertex(1).norm()) < eps)
                    face->line(line_no)->set_boundary_id(0);
            }
          else if (face_center_radius > middle_radius) // outer radius set boundary 1
            {
              face->set_boundary_id(1);
              for (const unsigned int line_no : face->line_indices())
                if (face->line(line_no)->at_boundary())
                  if (std::fabs(face->line(line_no)->vertex(0).norm() - face->line(line_no)->vertex(1).norm()) < eps)
                    face->line(line_no)->set_boundary_id(1);
            }
          else
            AssertThrow(false, ExcInternalError());
        }
  }

  template void
  colorize_quarter_hyper_shell<2>(Triangulation<2> &, const Point<2> &, const double, const double);
}

#endif
