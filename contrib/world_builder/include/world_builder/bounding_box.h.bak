// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef WORLD_BUILDER_BOUNDING_BOX_H
#define WORLD_BUILDER_BOUNDING_BOX_H


#include "assert.h"
#include "point.h"

#include <limits>


namespace WorldBuilder
{


  /**
   * A class that represents a box of arbitrary dimension <tt>spacedim</tt> and
   * with sides parallel to the coordinate axes, that is, a region
   * @f[
   * [x_0^L, x_0^U] \times ... \times [x_{spacedim-1}^L, x_{spacedim-1}^U],
   * @f]
   * where $(x_0^L , ..., x_{spacedim-1}^L)$ and $(x_0^U , ..., x_{spacedim-1}^U)$
   * denote the two vertices (bottom left and top right) which are used to
   * represent the box. The quantities $x_k^L$ and $x_k^U$ denote the "lower"
   * and "upper" bounds of values that are within the box for each coordinate
   * direction $k$.
   *
   * Geometrically, a bounding box is thus:
   * - 1D: a segment (represented by its vertices in the proper order)
   * - 2D: a rectangle (represented by the vertices V at bottom left, top right)
   * @code
   * .--------V
   * |        |
   * V--------.
   * @endcode
   *
   * - 3D: a cuboid (in which case the two vertices V follow the convention and
   * are not owned by the same face)
   * @code
   *   .------V
   *  /      /|
   * .------. |
   * |      | /
   * |      |/
   * V------.
   * @endcode
   *
   * Bounding boxes are, for example, useful in parallel distributed meshes to
   * give a general description of the owners of each portion of the mesh. More
   * generally, bounding boxes are often used to roughly describe a region of
   * space in which an object is contained; if a candidate point is not within
   * the bounding box (a test that is cheap to execute), then it is not necessary
   * to perform an expensive test whether the candidate point is in fact inside
   * the object itself. Bounding boxes are therefore often used as a first,
   * cheap rejection test before more detailed checks. As such, bounding boxes
   * serve many of the same purposes as the
   * [convex hull](https://en.wikipedia.org/wiki/Convex_hull), for which it is
   * also relatively straightforward to compute whether a point is inside or
   * outside, though not quite as cheap as for the bounding box.
   *
   * Taking the cross section of a BoundingBox<spacedim> orthogonal to a given
   * direction gives a box in one dimension lower: BoundingBox<spacedim - 1>.
   * In 3D, the 2 coordinates of the cross section of BoundingBox<3> can be
   * ordered in 2 different ways. That is, if we take the cross section orthogonal
   * to the y direction we could either order a 3D-coordinate into a
   * 2D-coordinate as $(x,z)$ or as $(z,x)$. This class uses the second
   * convention, corresponding to the coordinates being ordered cyclicly
   * $x \rightarrow y \rightarrow z \rightarrow x \rightarrow ... $
   * To be precise, if we take a cross section:
   *
   * | Orthogonal to | Cross section coordinates ordered as |
   * |:-------------:|:------------------------------------:|
   * |      x        |               (y, z)                 |
   * |      y        |               (z, x)                 |
   * |      z        |               (x, y)                 |
   *
   * This is according to the convention set by the function
   * <code>coordinate_to_one_dim_higher</code>.
   *
   * @note The majority of this class is copied from the deal.II library.
   */
  template <unsigned int spacedim>
  class BoundingBox
  {
    public:
      /**
       * Standard constructor. Creates an object that corresponds to a
       * box that corresponds to the entire space, i.e. a degenerate
       * box with end points at minus and plus infinity.
       */
      BoundingBox();

      /**
       * Standard constructor for non-empty boxes: it uses a pair of points
       * which describe the box: one for the bottom and one for the top
       * corner.
       */
      BoundingBox(const std::pair<Point<spacedim>, Point<spacedim>>
                  &boundary_points);

      /**
       * Construct the bounding box that encloses all the points in the given
       * container.
       *
       * The constructor supports any Container that provides begin() and end()
       * iterators to Point<spacedim> elements.
       */
      template <class Container>
      BoundingBox(const Container &points);

      /**
       * Return a reference to the boundary_points
       */
      std::pair<Point<spacedim>, Point<spacedim>> &get_boundary_points();

      /**
       * Return a const reference to the boundary_points
       */
      const std::pair<Point<spacedim>, Point<spacedim>> &get_boundary_points() const;

      /**
       * Test for equality.
       */
      bool
      operator==(const BoundingBox<spacedim> &box) const;

      /**
       * Test for inequality.
       */
      bool
      operator!=(const BoundingBox<spacedim> &box) const;

      /**
       * Enlarge the current object so that it contains @p other_bbox .
       * If the current object already contains @p other_bbox then it is not changed
       * by this function.
       */
      void
      merge_with(const BoundingBox<spacedim> &other_bbox);

      /**
       * Return true if the point is inside the Bounding Box, false otherwise. The
       * parameter @p tolerance is a factor by which the bounding box is enlarged
       * relative to the dimensions of the bounding box in order to determine in a
       * numerically robust way whether the point is inside.
       *
       * This function is a wrapper for the function point_inside_implementation
       * to also test point+2pi is the longtiude is smaller then zero and point-2pi
       * if the longitude is larger than zero.
       */
      bool
      point_inside(
        const Point<spacedim> &p,
        const double tolerance = std::numeric_limits<double>::epsilon()) const;

      /**
       * Increase (or decrease) the size of the bounding box by the given amount.
       * After calling this method, the lower left corner of the bounding box will
       * have each coordinate decreased by @p amount, and the upper right corner
       * of the bounding box will have each coordinate increased by @p amount.
       *
       * If you call this method with a negative number, and one of the axes of the
       * original bounding box is smaller than amount/2, the method will trigger
       * an assertion.
       */
      void
      extend(const double amount);

      /**
       * Compute the volume (i.e. the dim-dimensional measure) of the BoundingBox.
       */
      double
      volume() const;

      /**
       * Returns the point in the center of the box.
       */
      Point<spacedim>
      center() const;

      /**
       * Returns the side length of the box in @p direction.
       */
      double
      side_length(const unsigned int direction) const;

      /**
       * Return the lower bound of the box in @p direction.
       */
      double
      lower_bound(const unsigned int direction) const;

      /**
       * Return the upper bound of the box in @p direction.
       */
      double
      upper_bound(const unsigned int direction) const;

    private:
      /**
       * Return true if the point is inside the Bounding Box, false otherwise. The
       * parameter @p tolerance is a factor by which the bounding box is enlarged
       * relative to the dimensions of the bounding box in order to determine in a
       * numerically robust way whether the point is inside.
       *
       * This function is supposed to be only used by the point_inside function.
       */
      bool
      point_inside_implementation(
        const Point<spacedim> &p,
        const double tolerance = std::numeric_limits<double>::epsilon()) const;

      std::pair<Point<spacedim>, Point<spacedim>> boundary_points;

  };


  /*------------------------ Inline functions: BoundingBox --------------------*/

#ifndef DOXYGEN


  template <>
  inline BoundingBox<2>::BoundingBox()
    :
    boundary_points ({{
      -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max(),
      cartesian
    },
    {
      +std::numeric_limits<double>::max(),
      +std::numeric_limits<double>::max(),
      cartesian
    }
  })
  {}

  template <>
  inline BoundingBox<3>::BoundingBox()
    :
    boundary_points ({{
      -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max(),
      cartesian
    },
    {
      +std::numeric_limits<double>::max(),
      +std::numeric_limits<double>::max(),
      +std::numeric_limits<double>::max(),
      cartesian
    }
  })
  {}



  template <unsigned int spacedim>
  inline BoundingBox<spacedim>::BoundingBox(
    const std::pair<Point<spacedim>, Point<spacedim>>
    &boundary_points_)
    :
    boundary_points (boundary_points_)
  {
    // We check the Bounding Box is not degenerate
    for (unsigned int i = 0; i < spacedim; ++i)
      WBAssert(boundary_points.first[i] <= boundary_points.second[i],
               "Bounding Box can't be created: the points' "
               "order should be bottom left, top right!");
  }



  template <unsigned int spacedim>
  template <class Container>
  inline BoundingBox<spacedim>::BoundingBox(const Container &points)
  {
    // Use the default constructor in case points is empty instead of setting
    // things to +oo and -oo
    if (points.size() > 0)
      {
        auto &min = boundary_points.first;
        auto &max = boundary_points.second;
        std::fill(min.begin_raw(),
                  min.end_raw(),
                  std::numeric_limits<double>::infinity());
        std::fill(max.begin_raw(),
                  max.end_raw(),
                  -std::numeric_limits<double>::infinity());

        for (const Point<spacedim> &point : points)
          for (unsigned int d = 0; d < spacedim; ++d)
            {
              min[d] = std::min(min[d], point[d]);
              max[d] = std::max(max[d], point[d]);
            }
      }
  }



  template <unsigned int spacedim>
  inline std::pair<Point<spacedim>, Point<spacedim>> &BoundingBox<spacedim>::get_boundary_points()
  {
    return this->boundary_points;
  }



  template <unsigned int spacedim>
  inline const std::pair<Point<spacedim>, Point<spacedim>> &BoundingBox<spacedim>::get_boundary_points() const
  {
    return this->boundary_points;
  }



  template <unsigned int spacedim>
  inline bool
  BoundingBox<spacedim>::
  operator==(const BoundingBox<spacedim> &box) const
  {
    return boundary_points == box.boundary_points;
  }



  template <unsigned int spacedim>
  inline bool
  BoundingBox<spacedim>::
  operator!=(const BoundingBox<spacedim> &box) const
  {
    return boundary_points != box.boundary_points;
  }



  template <unsigned int spacedim>
  inline void
  BoundingBox<spacedim>::extend(const double amount)
  {
    for (unsigned int d = 0; d < spacedim; ++d)
      {
        boundary_points.first[d] -= amount;
        boundary_points.second[d] += amount;
        WBAssert(boundary_points.first[d] <= boundary_points.second[d],
                 "Bounding Box can't be shrunk this much: the points' "
                 "order should remain bottom left, top right.");
      }
  }


  template <unsigned int spacedim>
  inline
  bool
  BoundingBox<spacedim>::point_inside(const Point<spacedim> &point,
                                      const double tolerance) const
  {
    WBAssert(boundary_points.first.get_coordinate_system() == point.get_coordinate_system(),
             "Cannot compare two points which represent different coordinate systems.");
    WBAssert(boundary_points.second.get_coordinate_system() == point.get_coordinate_system(),
             "Cannot compare two points which represent different coordinate systems.");

    if (point.get_coordinate_system() == CoordinateSystem::spherical)
      {
        Point<spacedim> other_point = point;
        if (spacedim == 2)
          {
            other_point[0] += point[0] < 0 ? 2.0 * Consts::PI : -2.0 * Consts::PI;
          }
        else
          {
            // spacedim == 3 (rad,long,lat)
            other_point[1] += point[1] < 0 ? 2.0 * Consts::PI : -2.0 * Consts::PI;
          }

        return (point_inside_implementation(point, tolerance) ||
                point_inside_implementation(other_point, tolerance));
      }

    return point_inside_implementation(point, tolerance);

  }

  template <unsigned int spacedim>
  inline
  bool
  BoundingBox<spacedim>::point_inside_implementation(const Point<spacedim> &p,
                                                     const double tolerance) const
  {
    WBAssert(boundary_points.first.get_coordinate_system() == p.get_coordinate_system(),
             "Cannot compare two points which represent different coordinate systems.");
    WBAssert(boundary_points.second.get_coordinate_system() == p.get_coordinate_system(),
             "Cannot compare two points which represent different coordinate systems.");

    for (unsigned int i = 0; i < spacedim; ++i)
      {
        // Bottom left-top right convention: the point is outside if it's smaller
        // than the first or bigger than the second boundary point The bounding
        // box is defined as a closed set
        if ((p[i] < this->boundary_points.first[i] -
             tolerance * std::abs(this->boundary_points.second[i] -
                                  this->boundary_points.first[i])) ||
            (p[i] > this->boundary_points.second[i] +
             tolerance * std::abs(this->boundary_points.second[i] -
                                  this->boundary_points.first[i])))
          return false;
      }
    return true;
  }



  template <unsigned int spacedim>
  inline
  void
  BoundingBox<spacedim>::merge_with(
    const BoundingBox<spacedim> &other_bbox)
  {
    for (unsigned int i = 0; i < spacedim; ++i)
      {
        this->boundary_points.first[i] =
          std::min(this->boundary_points.first[i],
                   other_bbox.boundary_points.first[i]);
        this->boundary_points.second[i] =
          std::max(this->boundary_points.second[i],
                   other_bbox.boundary_points.second[i]);
      }
  }



  template <unsigned int spacedim>
  double
  BoundingBox<spacedim>::volume() const
  {
    double vol = 1.0;
    for (unsigned int i = 0; i < spacedim; ++i)
      vol *= (this->boundary_points.second[i] - this->boundary_points.first[i]);
    return vol;
  }



  template <unsigned int spacedim>
  inline
  double
  BoundingBox<spacedim>::lower_bound(const unsigned int direction) const
  {
    WBAssert(direction, spacedim);

    return boundary_points.first[direction];
  }



  template <unsigned int spacedim>
  inline
  double
  BoundingBox<spacedim>::upper_bound(const unsigned int direction) const
  {
    WBAssert(direction, spacedim);

    return boundary_points.second[direction];
  }



  template <unsigned int spacedim>
  Point<spacedim>
  BoundingBox<spacedim>::center() const
  {
    Point<spacedim> point = boundary_points.first; // initialize to inherit coordinate system
    for (unsigned int i = 0; i < spacedim; ++i)
      point[i] = .5 * (boundary_points.first[i] + boundary_points.second[i]);

    return point;
  }



  template <unsigned int spacedim>
  inline
  double
  BoundingBox<spacedim>::side_length(const unsigned int direction) const
  {
    WBAssert(direction < spacedim, "Invalid index");

    return boundary_points.second[direction] - boundary_points.first[direction];
  }


#endif // DOXYGEN

}


#endif
