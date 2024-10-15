/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_compat_h
#define _aspect_compat_h

#include <aspect/global.h>
// C++11 related includes.
#include <array>
#include <functional>
#include <memory>

namespace big_mpi
{

  using dealii::Utilities::MPI::broadcast;

}

// deal.II 9.6 introduces the new MGTransferMF class as a replacement
// for MGTransferMatrixFree. Instead of putting an ifdef in every place,
// do this in one central location:
#if !DEAL_II_VERSION_GTE(9,6,0)
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
namespace dealii
{
  template <int dim, class NumberType>
  using MGTransferMF = MGTransferMatrixFree<dim,NumberType>;
}
#endif



// deal.II versions up to 9.5 had a poorly designed interface of the
// SphericalManifold class that made it impossible for us to use.
// This file thus contains a copy of it.
#if !DEAL_II_VERSION_GTE(9,6,0)

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

namespace aspect
{
  using namespace dealii;

  /**
   * The deal.II class SphericalManifold has a design flaw that made it
   * impossible to derive from the class. This is fixed post-9.5,
   * see https://github.com/dealii/dealii/pull/16242 and
   * https://github.com/dealii/dealii/pull/16248, but we can't
   * use deal.II 9.5 and earlier for this class. The current class
   * here is therefore a copy of the fixed class.
   */
  template <int dim, int spacedim = dim>
  class SphericalManifold : public Manifold<dim, spacedim>
  {
    public:
      /**
       * Constructor.
       *
       * @param[in] center The center of the coordinate system. Defaults to the
       * origin.
       */
      SphericalManifold(const Point<spacedim> center = Point<spacedim>());

      /**
       * Make a clone of this Manifold object.
       */
      virtual std::unique_ptr<Manifold<dim, spacedim>>
      clone() const override;

      /**
       * Given any two points in space, first project them on the surface
       * of a sphere with unit radius, then connect them with a geodesic
       * and find the intermediate point, and finally rescale the final
       * radius so that the resulting one is the convex combination of the
       * starting radii.
       */
      virtual Point<spacedim>
      get_intermediate_point(const Point<spacedim> &p1,
                             const Point<spacedim> &p2,
                             const double           w) const override;

      /**
       * Compute the derivative of the get_intermediate_point() function
       * with parameter w equal to zero.
       */
      virtual Tensor<1, spacedim>
      get_tangent_vector(const Point<spacedim> &x1,
                         const Point<spacedim> &x2) const override;

      /**
       * @copydoc Manifold::normal_vector()
       */
      virtual Tensor<1, spacedim>
      normal_vector(
        const typename Triangulation<dim, spacedim>::face_iterator &face,
        const Point<spacedim> &p) const override;

      /**
       * Compute the normal vectors to the boundary at each vertex.
       */
      virtual void
      get_normals_at_vertices(
        const typename Triangulation<dim, spacedim>::face_iterator &face,
        typename Manifold<dim, spacedim>::FaceVertexNormals &face_vertex_normals)
      const override;

      /**
       * Compute a new set of points that interpolate between the given points @p
       * surrounding_points. @p weights is a table with as many columns as @p
       * surrounding_points.size(). The number of rows in @p weights must match
       * the length of @p new_points.
       *
       * This function is optimized to perform on a collection
       * of new points, by collecting operations that are not dependent on the
       * weights outside of the loop over all new points.
       *
       * The implementation does not allow for @p surrounding_points and
       * @p new_points to point to the same array, so make sure to pass different
       * objects into the function.
       */
      virtual void
      get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                     const Table<2, double>                 &weights,
                     ArrayView<Point<spacedim>> new_points) const override;

      /**
       * Return a point on the spherical manifold which is intermediate
       * with respect to the surrounding points.
       */
      virtual Point<spacedim>
      get_new_point(const ArrayView<const Point<spacedim>> &vertices,
                    const ArrayView<const double>          &weights) const override;

      /**
       * The center of the spherical coordinate system.
       */
      const Point<spacedim> center;

    private:
      /**
       * Return a point on the spherical manifold which is intermediate
       * with respect to the surrounding points. This function uses a linear
       * average of the directions to find an estimated point. It returns a pair
       * of radius and direction from the center point to the candidate point.
       */
      std::pair<double, Tensor<1, spacedim>>
      guess_new_point(const ArrayView<const Tensor<1, spacedim>> &directions,
                      const ArrayView<const double>              &distances,
                      const ArrayView<const double>              &weights) const;

      /**
       * This function provides an internal implementation of the get_new_points()
       * interface.
       *
       * It computes a new set of points that interpolate between the given points
       * @p
       * surrounding_points. @p weights is an array view with as many entries as @p
       * surrounding_points.size() times @p new_points.size().
       *
       * This function is optimized to perform on a collection
       * of new points, by collecting operations that are not dependent on the
       * weights outside of the loop over all new points.
       *
       * The implementation does not allow for @p surrounding_points and
       * @p new_points to point to the same array, so make sure to pass different
       * objects into the function.
       */
      void
      do_get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                        const ArrayView<const double>          &weights,
                        ArrayView<Point<spacedim>>              new_points) const;

      /**
       * A manifold description to be used for get_new_point in 2d.
       */
      const PolarManifold<spacedim> polar_manifold;
  };
}

#else

// For sufficiently new deal.II versions, we can use the deal.II class, but to
// avoid name clashes, we have to import the class into namespace aspect. Once
// we rely on these sufficiently new versions of deal.II, we can not only remove
// the code above, but also the following lines, and in all places where we
// reference 'aspect::SphericalManifold' simply use 'SphericalManifold' instead
// (which then refers to the deal.II class).

#include <deal.II/grid/manifold_lib.h>
namespace aspect
{
  using dealii::SphericalManifold;
}

#endif

#endif
