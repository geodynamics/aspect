/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__compat_h
#define __aspect__compat_h

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

/*
 * cross_product
 */
#if !DEAL_II_VERSION_GTE(8,4,0)
template <int dim>
dealii::Tensor<1,dim> cross_product_3d(const dealii::Tensor<1,dim> &a, const dealii::Tensor<1,dim> &b)
{
  Assert (dim==3, dealii::ExcInternalError());
  dealii::Tensor<1,dim> result;
  dealii::cross_product(result, a, b);
  return result;
}

template <int dim>
dealii::Tensor<1,dim> cross_product_2d(const dealii::Tensor<1,dim> &a)
{
  Assert (dim==2, dealii::ExcInternalError());
  dealii::Tensor<1,dim> result;
  dealii::cross_product(result, a);
  return result;
}
#endif


/*
 * MPI::min() functions
 */
#if !DEAL_II_VERSION_GTE(8,3,0)
namespace dealii
{
  namespace Utilities
  {
    namespace MPI
    {
      template <typename T>
      inline
      T min (const T &t,
             const MPI_Comm &mpi_communicator)
      {
        return -max(-t, mpi_communicator);
      }

      template <typename T>
      inline
      void min (const std::vector<T> &values,
                const MPI_Comm       &mpi_communicator,
                std::vector<T>       &minima)
      {
        minima.resize(values.size());
        for (unsigned int i=0; i<values.size(); ++i)
          minima[i] = -values[i];
        max(minima, mpi_communicator, minima);
        for (unsigned int i=0; i<values.size(); ++i)
          minima[i] = -minima[i];
      }
    }
  }
}
#endif

/*
 * Unique_ptr functionality that was introduced in deal.II 8.3 and replaces
 * auto_ptr in CXX17. We need this to silence deprecation warnings in new
 * compilers.
 */
#if DEAL_II_VERSION_GTE(8,3,0)
#include <deal.II/base/std_cxx11/unique_ptr.h>
#else
#ifdef DEAL_II_WITH_CXX11

#  include <memory>
namespace dealii
{
  namespace std_cxx11
  {
    using std::unique_ptr;
  }
}

#else

#include <boost/scoped_ptr.hpp>

namespace dealii
{
  namespace std_cxx11
  {
    /**
     * Implementation of a basic replacement for C++11's std::unique_ptr class.
     *
     * BOOST does not have a replacement for std::unique_ptr (because unique_ptr
     * requires move semantics that aren't available unless you have a C++11
     * compiler -- in which case you also have std::unique_ptr; see for example
     * http://stackoverflow.com/questions/2953530/unique-ptr-boost-equivalent)
     *
     * Consequently, we emulate the class by just wrapping a boost::scoped_ptr
     * in the cheapest possible way -- by just deriving from it and repeating
     * the basic constructors. Everything else is inherited from the scoped_ptr
     * class.
     *
     * There is no overhead to this approach: scoped_ptr cannot be copied or
     * moved. Instances of unique_ptr cannot be copied, and if you do not have a
     * C++11 compiler, then you cannot move anything anyway.
     */
    template <typename T>
    class unique_ptr : public boost::scoped_ptr<T>
    {
      public:
        unique_ptr () {}

        template<class Y>
        explicit unique_ptr (Y *p)
          :
          boost::scoped_ptr<T>(p)
        {}
    };

  }
}

#endif
#endif

/*
 * parallel::distributed::Triangulation::ghost_owners() function
 */
#if !DEAL_II_VERSION_GTE(8,4,0)

#include <deal.II/distributed/tria.h>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    template <int dim>
    std::set<types::subdomain_id>
    ghost_owners(const parallel::distributed::Triangulation<dim> &triangulation)
    {
      std::set<types::subdomain_id> neighbors;

      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell)
        {
          if (cell->is_ghost())
            {
              neighbors.insert(cell->subdomain_id());
            }
        }
      return neighbors;
    }
  }
}
#endif

#endif
