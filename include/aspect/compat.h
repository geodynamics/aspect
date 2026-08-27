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

#include <deal.II/base/config.h>
#include <deal.II/base/mpi.h>

// C++11 related includes.
#include <functional>
#include <memory>

#include <deal.II/grid/manifold_lib.h>

#if !DEAL_II_VERSION_GTE(9,7,0)
#  include <deal.II/grid/grid_generator.h>
#endif

#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

namespace aspect
{
  namespace big_mpi
  {

    using dealii::Utilities::MPI::broadcast;

  }

  // deal.II 9.6 introduced MGTransferMF as a replacement for
  // MGTransferMatrixFree; deal.II 9.9 renamed it back.
#if DEAL_II_VERSION_GTE(9,9,0)
  template <int dim, typename NumberType>
  using MGTransferType = dealii::MGTransferMatrixFree<dim, NumberType>;
#else
  template <int dim, typename NumberType>
  using MGTransferType = dealii::MGTransferMF<dim, NumberType>;
#endif

  // Global-coarsening transfers share the same class as local smoothing
  // from deal.II 9.7 onward; 9.6 still needs MGTransferGlobalCoarsening.
#if DEAL_II_VERSION_GTE(9,7,0)
  template <int dim, typename NumberType>
  using GCMGTransferType = MGTransferType<dim, NumberType>;
#else
  template <int dim, typename NumberType>
  using GCMGTransferType = dealii::MGTransferGlobalCoarsening<dim, dealii::LinearAlgebra::distributed::Vector<NumberType>>;
#endif


#if !DEAL_II_VERSION_GTE(9,7,0)
  /**
   * A type alias for the SmartPointer class that makes sure the new
   * name of the class, ObserverPointer, can be used in all versions of
   * deal.II. The SmartPointer class was renamed to ObserverPointer in
   * deal.II 9.7.
   */
  template <typename T, typename P = void>
  using ObserverPointer = dealii::SmartPointer<T, P>;

  /**
   * Same for the transition from Subscriptor to EnableObserverPointer.
   */
  using EnableObserverPointer = dealii::Subscriptor;
#endif


  using dealii::SphericalManifold;



// deal.II version 9.7 introduces a new class VectorFunctionFromTensorFunctionObject
// that we would like to use also for earlier versions
#if !DEAL_II_VERSION_GTE(9,7,0)

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

#endif

// deal.II versions up to 9.6 had a bug for very thin shell geometries.
// This function contains a fixed version.
#if !DEAL_II_VERSION_GTE(9,7,0)

  template <int dim>
  void
  colorize_quarter_hyper_shell(Triangulation<dim>  &tria,
                               const Point<dim>   &center,
                               const double      inner_radius,
                               const double      outer_radius);
#endif

  // deal.II 9.8 made ReferenceCell a template class, whereas older versions
  // had it as a non-template class. This is a problem.
  // Rather than litter our own code base with #ifdefs, we can just define the
  // templated class variant here for older deal.II versions, and then we can
  // use the same code in all versions.
#if !DEAL_II_VERSION_GTE(9,8,0)
  template <int dim> using ReferenceCell = dealii::ReferenceCell;
#endif

}

#endif
