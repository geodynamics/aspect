/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#ifndef _aspect_compat_h
#define _aspect_compat_h

#include <aspect/global.h>
/*
 * Fixed colorization of parallelepiped
 */
#if !DEAL_II_VERSION_GTE(8,4,0)
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
namespace dealii
{
  namespace GridGenerator
  {
    // Parallelepiped implementation in 1d, 2d, and 3d. @note The
    // implementation in 1d is similar to hyper_rectangle(), and in 2d is
    // similar to parallelogram().
    //
    // The GridReordering::reorder_grid is made use of towards the end of
    // this function. Thus the triangulation is explicitly constructed for
    // all dim here since it is slightly different in that respect
    // (cf. hyper_rectangle(), parallelogram()).
    template<int dim>
    void
    subdivided_parallelepiped (Triangulation<dim>  &tria,
#ifndef _MSC_VER
                               const unsigned int(&n_subdivisions)[dim],
#else
                               const unsigned int *n_subdivisions,
#endif
                               const Point<dim>   (&corners) [dim],
                               const bool           colorize)
    {
      // Zero n_subdivisions is the origin only, which makes no sense
      for (unsigned int i=0; i<dim; ++i)
        Assert (n_subdivisions[i]>0, ExcInvalidRepetitions(n_subdivisions[i]));

      // Check corners do not overlap (unique)
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=i+1; j<dim; ++j)
          Assert ((corners[i]!=corners[j]),
                  ExcMessage ("Invalid distance between corner points of parallelepiped."));

      // Create a list of points
      Point<dim> delta[dim];

      for (unsigned int i=0; i<dim; ++i)
        delta[i] = corners[i]/n_subdivisions[i];
      std::vector<Point<dim> > points;

      switch (dim)
        {
          case 1:
            for (unsigned int x=0; x<=n_subdivisions[0]; ++x)
              points.push_back (double(x)*delta[0]);
            break;

          case 2:
            for (unsigned int y=0; y<=n_subdivisions[1]; ++y)
              for (unsigned int x=0; x<=n_subdivisions[0]; ++x)
                points.push_back (double(x)*delta[0] + double(y)*delta[1]);
            break;

          case 3:
            for (unsigned int z=0; z<=n_subdivisions[2]; ++z)
              for (unsigned int y=0; y<=n_subdivisions[1]; ++y)
                for (unsigned int x=0; x<=n_subdivisions[0]; ++x)
                  points.push_back (double(x)*delta[0] + double(y)*delta[1] + double(z)*delta[2]);
            break;

          default:
            Assert (false, ExcNotImplemented());
        }

      // Prepare cell data
      unsigned int n_cells = 1;
      for (unsigned int i=0; i<dim; ++i)
        n_cells *= n_subdivisions[i];
      std::vector<CellData<dim> > cells (n_cells);

      // Create fixed ordering of
      switch (dim)
        {
          case 1:
            for (unsigned int x=0; x<n_subdivisions[0]; ++x)
              {
                cells[x].vertices[0] = x;
                cells[x].vertices[1] = x+1;

                // wipe material id
                cells[x].material_id = 0;
              }
            break;

          case 2:
          {
            // Shorthand
            const unsigned int n_dy = n_subdivisions[1];
            const unsigned int n_dx = n_subdivisions[0];

            for (unsigned int y=0; y<n_dy; ++y)
              for (unsigned int x=0; x<n_dx; ++x)
                {
                  const unsigned int    c = y*n_dx         + x;
                  cells[c].vertices[0] = y*(n_dx+1)     + x;
                  cells[c].vertices[1] = y*(n_dx+1)     + x+1;
                  cells[c].vertices[2] = (y+1)*(n_dx+1) + x;
                  cells[c].vertices[3] = (y+1)*(n_dx+1) + x+1;

                  // wipe material id
                  cells[c].material_id = 0;
                }
          }
          break;

          case 3:
          {
            // Shorthand
            const unsigned int n_dz = n_subdivisions[2];
            const unsigned int n_dy = n_subdivisions[1];
            const unsigned int n_dx = n_subdivisions[0];

            for (unsigned int z=0; z<n_dz; ++z)
              for (unsigned int y=0; y<n_dy; ++y)
                for (unsigned int x=0; x<n_dx; ++x)
                  {
                    const unsigned int    c = z*n_dy*n_dx             + y*n_dx         + x;

                    cells[c].vertices[0] = z*(n_dy+1)*(n_dx+1)     + y*(n_dx+1)     + x;
                    cells[c].vertices[1] = z*(n_dy+1)*(n_dx+1)     + y*(n_dx+1)     + x+1;
                    cells[c].vertices[2] = z*(n_dy+1)*(n_dx+1)     + (y+1)*(n_dx+1) + x;
                    cells[c].vertices[3] = z*(n_dy+1)*(n_dx+1)     + (y+1)*(n_dx+1) + x+1;
                    cells[c].vertices[4] = (z+1)*(n_dy+1)*(n_dx+1) + y*(n_dx+1)     + x;
                    cells[c].vertices[5] = (z+1)*(n_dy+1)*(n_dx+1) + y*(n_dx+1)     + x+1;
                    cells[c].vertices[6] = (z+1)*(n_dy+1)*(n_dx+1) + (y+1)*(n_dx+1) + x;
                    cells[c].vertices[7] = (z+1)*(n_dy+1)*(n_dx+1) + (y+1)*(n_dx+1) + x+1;

                    // wipe material id
                    cells[c].material_id = 0;
                  }
            break;
          }

          default:
            Assert (false, ExcNotImplemented());
        }

      // Create triangulation
      // reorder the cells to ensure that they satisfy the convention for
      // edge and face directions
      GridReordering<dim>::reorder_cells(cells, true);
      tria.create_triangulation (points, cells, SubCellData());

      // Finally assign boundary indicators according to hyper_rectangle
      if (colorize)
        {
          typename Triangulation<dim>::active_cell_iterator
          cell = tria.begin_active(),
          endc = tria.end();
          for (; cell!=endc; ++cell)
            {
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                {
                  if (cell->face(face)->at_boundary())
#if DEAL_II_VERSION_GTE(8,3,0)
                    cell->face(face)->set_boundary_id(face);
#else
                    cell->face(face)->set_boundary_indicator (face);
#endif
                }
            }
        }
    }
  }
}
#endif

// C++11 related includes. Can be removed when we require C++11.
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/std_cxx11/bind.h>
#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>


// We would like to use a function from SolverControl that was introduced after
// deal.II 8.5. For older versions use this derived class instead that implements
// the function and the typedef guarantees it is found before deal.II's class.
#if !DEAL_II_VERSION_GTE(9,0,0)

#include <deal.II/lac/solver_control.h>

namespace aspect
{
  using namespace dealii;

  class SolverControl : public dealii::SolverControl
  {
    public:
      SolverControl(const unsigned int n           = 100,
                    const double       tol         = 1.e-10,
                    const bool         log_history = false,
                    const bool         log_result  = true)
        :
        dealii::SolverControl (n, tol, log_history, log_result)
      {}

      dealii::SolverControl::State
      check (const unsigned int step,
             const double check_value)
      {
        dealii::SolverControl::State return_value = dealii::SolverControl::check(step, check_value);

        if (step == 0)
          history_data.resize(history_data.size()+1);
        return return_value;
      }


      const std::vector<double> &get_history_data() const
      {
        Assert (history_data_enabled, ExcHistoryDataRequired());
        Assert (history_data.size() > 0,
                ExcMessage("The SolverControl object was asked for the solver history "
                           "data, but there is no data. Possibly you requested the data before the "
                           "solver was run."));

        return history_data;
      }
  };
}
#endif


#endif
