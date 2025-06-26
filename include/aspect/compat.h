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

namespace aspect
{
  namespace big_mpi
  {

    using dealii::Utilities::MPI::broadcast;

  }
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


// deal.II version 9.7 introduces a new class VectorFunctionFromTensorFunctionObject
// that we would like to use also for earlier versions
#if !DEAL_II_VERSION_GTE(9,7,0)

namespace aspect
{
  using namespace dealii;

  /**
   * This class is built as a means of translating the <code>Tensor<1,dim,
   * RangeNumberType> </code> values produced by function objects that
   * for a given point return a tensor,
   * and returning them as a multiple component version of the same thing as a
   * Vector for use in, for example, the VectorTools::interpolate or the many
   * other functions taking Function objects. It allows the user to place the
   * desired components into an <tt>n_components</tt> long vector starting at
   * the <tt>selected_component</tt> location in that vector and have all other
   * components be 0.
   *
   * For example: Say you created a function object that returns the gravity
   * (here, a radially inward pointing vector of magnitude 9.81):
   * @code
   *   const auto gravity
   *     = [](const Point<dim> &p) -> Tensor<1,dim> { return -9.81 * (p /
   * p.norm()); }
   * @endcode
   * and you want to interpolate this onto your mesh using the
   * VectorTools::interpolate function, with a finite element for the
   * DoFHandler object has 3 copies of a finite element with <tt>dim</tt>
   * components, for a total of 3*dim components. To interpolate onto that
   * DoFHandler, you need an object of type Function that has 3*dim vector
   * components. Creating such an object from the existing <code>gravity</code>
   * object is done using this piece of code:
   * @code
   * VectorFunctionFromTensorFunctionObject<dim, RangeNumberType>
   *   gravity_vector_function(gravity, 0, 3*dim);
   * @endcode
   * where the <code>dim</code> components of the `gravity` function are placed
   * into the first <code>dim</code> components of the function object.
   *
   * @ingroup functions
   */
  template <int dim, typename RangeNumberType = double>
  class VectorFunctionFromTensorFunctionObject
    : public Function<dim, RangeNumberType>
  {
    public:
      /**
       * Given a function object that takes a <tt>Point</tt> and returns a
       * <tt>Tensor<1,dim, RangeNumberType></tt> value, convert this into an object
       * that matches the Function@<dim@> interface.
       *
       * By default, create a Vector object of the same size as
       * <tt>tensor_function</tt> returns, i.e., with <tt>dim</tt> components.
       *
       * @param tensor_function_object The TensorFunction that will form `dim` components of
       * the resulting Vector Function object.
       * @param n_components The total number of vector components of the
       * resulting TensorFunction object. This clearly has to be at least `dim`.
       * @param selected_component The first component that should be filled by
       * the first argument.  This should be such that the entire tensor_function
       * fits inside the <tt>n_component</tt> length return vector.
       */
      explicit VectorFunctionFromTensorFunctionObject(
        const std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>
        &tensor_function_object,
        const unsigned int selected_component = 0,
        const unsigned int n_components       = dim);

      /**
       * This destructor is defined as virtual so as to coincide with all other
       * aspects of class.
       */
      virtual ~VectorFunctionFromTensorFunctionObject() override = default;

      /**
       * Return a single component of a vector-valued function at a given point.
       */
      virtual RangeNumberType
      value(const Point<dim> &p, const unsigned int component = 0) const override;

      /**
       * Return all components of a vector-valued function at a given point.
       *
       * <tt>values</tt> shall have the right size beforehand, i.e. #n_components.
       */
      virtual void
      vector_value(const Point<dim>        &p,
                   Vector<RangeNumberType> &values) const override;

      /**
       * Return all components of a vector-valued function at a list of points.
       *
       * <tt>value_list</tt> shall be the same size as <tt>points</tt> and each
       * element of the vector will be passed to vector_value() to evaluate the
       * function
       */
      virtual void
      vector_value_list(
        const std::vector<Point<dim>>        &points,
        std::vector<Vector<RangeNumberType>> &value_list) const override;

    private:
      /**
       * The TensorFunction object which we call when this class's vector_value()
       * or vector_value_list() functions are called.
       */
      const std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>
      tensor_function_object;

      /**
       * The first vector component whose value is to be filled by the given
       * TensorFunction.  The values will be placed in components
       * selected_component to selected_component+dim-1 for a
       * <tt>TensorFunction<1,dim, RangeNumberType></tt> object.
       */
      const unsigned int selected_component;
  };


  template <int dim, typename RangeNumberType>
  VectorFunctionFromTensorFunctionObject<dim, RangeNumberType>::
  VectorFunctionFromTensorFunctionObject(
    const std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>
    &tensor_function_object,
    const unsigned int selected_component,
    const unsigned int n_components)
    : Function<dim, RangeNumberType>(n_components)
    , tensor_function_object(tensor_function_object)
    , selected_component(selected_component)
  {
    // Verify that the Tensor<1,dim,RangeNumberType> will fit in the given length
    // selected_components and not hang over the end of the vector.
    AssertIndexRange(selected_component + dim - 1, this->n_components);
  }



  template <int dim, typename RangeNumberType>
  inline RangeNumberType
  VectorFunctionFromTensorFunctionObject<dim, RangeNumberType>::value(
    const Point<dim>  &p,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);

    // if the requested component is out of the range selected, then we can
    // return early
    if ((component < selected_component) ||
        (component >= selected_component + dim))
      return 0;

    // otherwise retrieve the values from the <tt>tensor_function</tt> to be
    // placed at the <tt>selected_component</tt> to
    // <tt>selected_component + dim - 1</tt> elements of the <tt>Vector</tt>
    // values and pick the correct one
    const Tensor<1, dim, RangeNumberType> tensor_value =
      tensor_function_object(p);

    return tensor_value[component - selected_component];
  }


  template <int dim, typename RangeNumberType>
  inline void
  VectorFunctionFromTensorFunctionObject<dim, RangeNumberType>::vector_value(
    const Point<dim>        &p,
    Vector<RangeNumberType> &values) const
  {
    Assert(values.size() == this->n_components,
           ExcDimensionMismatch(values.size(), this->n_components));

    // Retrieve the values from the <tt>tensor_function</tt> to be placed at
    // the <tt>selected_component</tt> to
    // <tt>selected_component + dim - 1</tt> elements of the <tt>Vector</tt>
    // values.
    const Tensor<1, dim, RangeNumberType> tensor_value =
      tensor_function_object(p);

    // First we make all elements of values = 0
    values = 0;

    // Second we adjust the desired components to take on the values in
    // <tt>tensor_value</tt>.
    for (unsigned int i = 0; i < dim; ++i)
      values(i + selected_component) = tensor_value[i];
  }


  /**
   * Member function <tt>vector_value_list </tt> is the interface for giving a
   * list of points (<code>vector<Point<dim>></code>) of which to evaluate
   * using the <tt>vector_value</tt> member function.  Again, this function is
   * written so as to not replicate the function definition but passes each
   * point on to <tt>vector_value</tt> to be evaluated.
   */
  template <int dim, typename RangeNumberType>
  void
  VectorFunctionFromTensorFunctionObject<dim, RangeNumberType>::vector_value_list(
    const std::vector<Point<dim>>        &points,
    std::vector<Vector<RangeNumberType>> &value_list) const
  {
    Assert(value_list.size() == points.size(),
           ExcDimensionMismatch(value_list.size(), points.size()));

    const unsigned int n_points = points.size();

    for (unsigned int p = 0; p < n_points; ++p)
      VectorFunctionFromTensorFunctionObject<dim, RangeNumberType>::vector_value(
        points[p], value_list[p]);
  }
}

#endif

// deal.II versions up to 9.6 had a bug for very thin shell geometries.
// This function contains a fixed version.
#if !DEAL_II_VERSION_GTE(9,7,0)

#include <deal.II/grid/grid_generator.h>

namespace aspect
{
  template <int dim>
  void
  colorize_quarter_hyper_shell(Triangulation<dim>  &tria,
                               const Point<dim>   &center,
                               const double      inner_radius,
                               const double      outer_radius);
}
#endif

#endif
