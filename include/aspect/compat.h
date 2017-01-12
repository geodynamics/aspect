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

//#include <aspect/global.h>

#include <deal.II/grid/tria_boundary_lib.h>

namespace dealii
{
  template <int dim>
  class SpecialConeBoundary : public StraightBoundary<dim>
  {
    public:

      /**
       * Constructor. Here the boundary object is constructed. The points
       * <tt>x_0</tt> and <tt>x_1</tt> describe the starting and ending points of
       * the axis of the (truncated) cone. <tt>radius_0</tt> denotes the radius
       * corresponding to <tt>x_0</tt> and <tt>radius_1</tt> the one corresponding
       * to <tt>x_1</tt>.
       */
      SpecialConeBoundary (const double radius_0,
                           const double radius_1,
                           const Point<dim> x_0,
                           const Point<dim> x_1);

      /**
       * Return the radius of the (truncated) cone at given point <tt>x</tt> on
       * the axis.
       */
      double get_radius (const Point<dim> x) const;

      /**
       * Refer to the general documentation of this class and the documentation of
       * the base class.
       */
      virtual
      Point<dim>
      get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

      /**
       * Refer to the general documentation of this class and the documentation of
       * the base class.
       */
      virtual
      Point<dim>
      get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;

      /**
       * Refer to the general documentation of this class and the documentation of
       * the base class.
       *
       * Calls @p get_intermediate_points_between_points.
       */
      virtual
      void
      get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                       std::vector<Point<dim> > &points) const;

      virtual
      Tensor<1,dim>
      normal_vector (const typename Triangulation<dim>::face_iterator &face,
                     const Point<dim> &p) const;


      /**
       * Refer to the general documentation of this class and the documentation of
       * the base class.
       *
       * Only implemented for <tt>dim=3</tt> and for <tt>points.size()==1</tt>.
       */
      virtual
      void
      get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
                                       std::vector<Point<dim> > &points) const;

      /**
       * Compute the normals to the boundary at the vertices of the given face.
       *
       * Refer to the general documentation of this class and the documentation of
       * the base class.
       */
      virtual
      void
      get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                               typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const;

    protected:
      /**
       * First radius of the (truncated) cone.
       */
      const double radius_0;

      /**
       * Second radius of the (truncated) cone.
       */
      const double radius_1;

      /**
       * Starting point of the axis.
       */
      const Point<dim> x_0;

      /**
       * Ending point of the axis.
       */
      const Point<dim> x_1;

    private:
      /**
       * Called by @p get_intermediate_points_on_line and by @p
       * get_intermediate_points_on_quad.
       *
       * Refer to the general documentation of @p get_intermediate_points_on_line
       * in the documentation of the base class.
       */
      void
      get_intermediate_points_between_points (const Point<dim> &p0,
                                              const Point<dim> &p1,
                                              std::vector<Point<dim> > &points) const;
  };




}
namespace dealii
{
  template<int dim>
  SpecialConeBoundary<dim>::SpecialConeBoundary (const double radius_0,
                                                 const double radius_1,
                                                 const Point<dim> x_0,
                                                 const Point<dim> x_1)
    :
    radius_0 (radius_0),
    radius_1 (radius_1),
    x_0 (x_0),
    x_1 (x_1)
  {}



  template<int dim>
  double SpecialConeBoundary<dim>::get_radius (Point<dim> x) const
  {
    for (unsigned int i = 0; i < dim; ++i)
      if ((x_1 (i) - x_0 (i)) != 0)
        return (radius_1 - radius_0) * (x (i) - x_0 (i)) / (x_1 (i) - x_0 (i)) + radius_0;

    return 0;
  }



  template<int dim>
  void
  SpecialConeBoundary<dim>::
  get_intermediate_points_between_points (const Point<dim> &p0,
                                          const Point<dim> &p1,
                                          std::vector<Point<dim> > &points) const
  {
    const unsigned int n = points.size ();
    const Tensor<1,dim> axis = x_1 - x_0;

    Assert (n > 0, ExcInternalError ());

    const std::vector<Point<1> > &line_points = this->get_line_support_points(n);

    for (unsigned int i=0; i<n; ++i)
      {
        const double x = line_points[i+1][0];

        // Compute the current point.
        const Point<dim> x_i = (1-x)*p0 + x*p1;
        // To project this point on the boundary of the cone we first compute
        // the orthogonal projection of this point onto the axis of the cone.
        const double c = (x_i - x_0) * axis / (axis*axis);
        const Point<dim> x_ip = x_0 + c * axis;
        // Compute the projection of the middle point on the boundary of the
        // cone.
        points[i] = x_ip + get_radius (x_ip) *  (x_i - x_ip) / (x_i - x_ip).norm ();
      }
  }

  template<int dim>
  Point<dim>
  SpecialConeBoundary<dim>::
  get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
  {
    const Tensor<1,dim> axis = x_1 - x_0;
    // Compute the middle point of the line.
    const Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);
    // To project it on the boundary of the cone we first compute the orthogonal
    // projection of the middle point onto the axis of the cone.
    const double c = (middle - x_0) * axis / (axis*axis);
    const Point<dim> middle_p = x_0 + c * axis;
    // Compute the projection of the middle point on the boundary of the cone.
    return middle_p + get_radius (middle_p) * (middle - middle_p) / (middle - middle_p).norm ();
  }



  template <>
  Point<3>
  SpecialConeBoundary<3>::
  get_new_point_on_quad (const Triangulation<3>::quad_iterator &quad) const
  {
    const int dim = 3;

    const Tensor<1,dim> axis = x_1 - x_0;
    // Compute the middle point of the quad.
    const Point<dim> middle = StraightBoundary<3,3>::get_new_point_on_quad (quad);
    // Same algorithm as above: To project it on the boundary of the cone we
    // first compute the orthogonal projection of the middle point onto the axis
    // of the cone.
    const double c = (middle - x_0) * axis / (axis*axis);
    const Point<dim> middle_p = x_0 + c * axis;
    // Compute the projection of the middle point on the boundary of the cone.
    return middle_p + get_radius (middle_p) * (middle - middle_p) / (middle - middle_p).norm ();
  }



//template<int dim>
//Point<dim>
//SpecialConeBoundary<dim>::
//get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &) const
//{
//  Assert (false, ExcImpossibleInDim (dim));
//
//  return Point<dim>();
//}
//


  template<int dim>
  void
  SpecialConeBoundary<dim>::
  get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                   std::vector<Point<dim> > &points) const
  {
    if (points.size () == 1)
      points[0] = get_new_point_on_line (line);
    else
      get_intermediate_points_between_points (line->vertex (0), line->vertex (1), points);
  }


  template<>
  void
  SpecialConeBoundary<3>::
  get_intermediate_points_on_quad (const Triangulation<3>::quad_iterator &quad,
                                   std::vector<Point<3> > &points) const
  {
    if (points.size () == 1)
      points[0] = get_new_point_on_quad (quad);
    else
      {
        unsigned int n = static_cast<unsigned int> (std::sqrt (static_cast<double> (points.size ())));

        Assert (points.size () == n * n, ExcInternalError ());

        std::vector<Point<3> > points_line_0 (n);
        std::vector<Point<3> > points_line_1 (n);

        get_intermediate_points_on_line (quad->line (0), points_line_0);
        get_intermediate_points_on_line (quad->line (1), points_line_1);

        std::vector<Point<3> > points_line_segment (n);

        for (unsigned int i = 0; i < n; ++i)
          {
            get_intermediate_points_between_points (points_line_0[i],
                                                    points_line_1[i],
                                                    points_line_segment);

            for (unsigned int j = 0; j < n; ++j)
              points[i * n + j] = points_line_segment[j];
          }
      }
  }



//template <int dim>
//void
//SpecialConeBoundary<dim>::
//get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &,
//                                 std::vector<Point<dim> > &) const
//{
//  Assert (false, ExcImpossibleInDim (dim));
//}
//
//
//template<>
//void
//SpecialConeBoundary<1>::
//get_normals_at_vertices (const Triangulation<1>::face_iterator &,
//                         Boundary<1,1>::FaceVertexNormals &) const
//{
//  Assert (false, ExcImpossibleInDim (1));
//}



  template<int dim>
  void
  SpecialConeBoundary<dim>::
  get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                           typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
  {
    const Tensor<1,dim> axis = x_1 - x_0;

    for (unsigned int vertex = 0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
      {
        // Compute the orthogonal projection of the vertex onto the axis of the
        // cone.
        const double c = (face->vertex (vertex) - x_0) * axis / (axis*axis);
        const Point<dim> vertex_p = x_0 + c * axis;
        // Then compute the vector pointing from the point <tt>vertex_p</tt> on
        // the axis to the vertex.
        const Tensor<1,dim> axis_to_vertex = face->vertex (vertex) - vertex_p;

        face_vertex_normals[vertex] = axis_to_vertex / axis_to_vertex.norm ();
      }
  }


  template<int dim>
  Tensor<1,dim>
  SpecialConeBoundary<dim>::
  normal_vector (const typename Triangulation<dim>::face_iterator &,
                 const Point<dim> &p) const
  {
    // TODO only for cone opening along z-axis
    AssertThrow (dim == 3, ExcInternalError());
    AssertThrow (radius_0 == 0., ExcInternalError());
    AssertThrow (x_0[0] == 0., ExcInternalError());
    const double c_squared = (radius_1 / x_1[dim-1])*(radius_1 / x_1[dim-1]);
    Tensor<1,dim> normal = p;
    normal[0] *= -2.0/c_squared;
    normal[1] *= -2.0/c_squared;
    normal[dim-1] *= 2.0;

    return normal/normal.norm();
  }
}

#endif
