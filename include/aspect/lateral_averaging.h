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


#ifndef _aspect_lateral_averaging_h
#define _aspect_lateral_averaging_h

#include <aspect/simulator_access.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  using namespace dealii;

  namespace internal
  {
    /**
     * This is the base class for all the functors implemented. Each is
     * used to compute one of the properties that will be laterally
     * averaged. The point of the base class is to allow handing over
     * a variable of type
     * <code> std::vector<std::unique_ptr<FunctorBase<dim> > > </code> to the
     * LateralAveraging::get_averages() function.
     */
    template <int dim>
    class FunctorBase
    {
      public:
        /**
         * operator() will have @p in and @p out filled out if @p true. By default
         * returns false.
         */
        virtual
        bool
        need_material_properties() const;

        /**
         * If this functor needs additional material model outputs create them
         * in here. By default, this does nothing.
         */
        virtual
        void
        create_additional_material_model_outputs (const unsigned int n_points,
                                                  MaterialModel::MaterialModelOutputs<dim> &outputs) const;

        /**
         * Called once at the beginning of compute_lateral_averages() to setup
         * internal data structures with the number of quadrature points.
         */
        virtual
        void
        setup(const unsigned int q_points);

        /**
         * This takes @p in material model inputs and @p out outputs (which are filled
         * if need_material_properties() == true), an initialized FEValues
         * object for a cell, and the current solution vector as inputs.
         * Functions in derived classes should then evaluate the desired quantity
         * and return the results in the output vector, which is q_points long.
         */
        virtual
        void
        operator()(const MaterialModel::MaterialModelInputs<dim> &in,
                   const MaterialModel::MaterialModelOutputs<dim> &out,
                   const FEValues<dim> &fe_values,
                   const LinearAlgebra::BlockVector &solution,
                   std::vector<double> &output) = 0;

        /**
         * Provide an (empty) virtual destructor.
         */
        virtual
        ~FunctorBase();
    };
  }

  /**
   * LateralAveraging is a class that performs various averaging operations
   * on the solution.  The functions of this class take a vector as an argument.
   * The model is divided into depth slices where the number of slices is the
   * length of the output vector. Each function averages a specific quantity
   * (as specified by their name), and that quantity is averaged laterally
   * for each depth slice.

   * Plugins may access the LateralAveraging plugin through the SimulatorAccess
   * function get_lateral_averaging(), and then query that for the desired
   * averaged quantity.
   *
   * @ingroup Simulator
   */
  template <int dim>
  class LateralAveraging : public SimulatorAccess<dim>
  {
    public:
      /**
       * @deprecated: This function is deprecated and only maintained for backward compatibilty.
       * Use the function compute_lateral_averages() with the same arguments instead.
       */
      DEAL_II_DEPRECATED
      std::vector<std::vector<double> >
      get_averages(const unsigned int n_slices,
                   const std::vector<std::string> &property_names) const;

      /**
       * Return a depth profile of lateral averages of the selected
       * @p property_names. This function is a convenience interface for
       * the other functions of the same name and is more efficient for
       * multiple properties than calling multiple of the single property functions
       * sequentially. This version of the function creates @p n_slices
       * equidistant depth slices to compute averages in.
       *
       * @param n_slices The number of equidistant depth slices
       * to perform the averaging in.
       * @param property_names Names of the available quantities to average.
       * Check the implementation of this function for available names.
       * @return The output vector of laterally averaged values. Each vector
       * has the same size of @p n_slices, and there are
       * as many vectors returned as names in @p property_names.
       */
      std::vector<std::vector<double> >
      compute_lateral_averages(const unsigned int n_slices,
                               const std::vector<std::string> &property_names) const;

      /**
       * Return a depth profile of lateral averages of the selected
       * @p property_names. This function is a convenience interface for
       * the other functions of the same name and is more efficient for
       * multiple properties than calling multiple of the single property
       * functions sequentially. This version of the function uses @p depth_bounds
       * to determine the upper and lower boundary of each depth slice,
       * and returns one averaged value per slice (= one less than the
       * number of depth bounds).
       *
       * @param depth_bounds The boundaries of the depth slices to compute
       * averages in. It is expected to be a vector of monotonically increasing
       * depth values with consecutive entries representing the minimum and
       * maximum depth of each depth slice. All depths smaller than
       * entry 0 or larger than the last entry are ignored, and depths
       * between entries 0 and 1 fall into depth slice 0, and so on.
       * @param property_names Names of the available quantities to average.
       * Check the implementation of this function for available names.
       * @return The output vector of laterally averaged values. The length
       * of each vector is one less than the number of @p depth_bounds,
       * and there are as many vectors returned as names in @p property_names.
       */
      std::vector<std::vector<double> >
      compute_lateral_averages(const std::vector<double> &depth_bounds,
                               const std::vector<std::string> &property_names) const;

      /**
       * Return a depth profile of lateral averages. This function is the
       * low-level implementation for the other functions of the same name
       * and is more efficient for multiple properties than calling multiple
       * of the single property functions sequentially. This version of the
       * function uses @p depth_bounds to determine the upper and lower
       * boundary of each depth slice, and returns one averaged value per slice
       * (= one less than the number of depth bounds).
       * The vector of functors @p functors must contain one or more
       * objects of classes that are derived from FunctorBase and are used to
       * fill the values vectors. All of the other functions in this class use
       * this low-level implementation for the actual computation.
       * Using this function allows to provide user-defined functors to average
       * properties not already implemented (e.g. to average additional
       * material model outputs).
       *
       * @param depth_bounds The boundaries of the depth slices to compute
       * averages in. It is expected to be a vector of monotonically increasing
       * depth values with consecutive entries representing the minimum and
       * maximum depth of each depth slice. All depths smaller than
       * entry 0 or larger than the last entry are ignored, and depths
       * between entries 0 and 1 fall into depth slice 0, and so on.
       * @param functors Instances of a class derived from FunctorBase
       * that are used to compute the averaged properties.
       * @return The output vectors of depth averaged values. The
       * function returns one vector of doubles per property and uses
       * @p depth_bounds to determine the depth extents of each slice.
       * Each returned vector has the same size (one entry less than
       * the number of @p depth_bounds).
       */
      std::vector<std::vector<double> >
      compute_lateral_averages(const std::vector<double> &depth_bounds,
                               std::vector<std::unique_ptr<internal::FunctorBase<dim> > > &functors) const;

      /**
       * Fill the argument with a set of lateral averages of the current
       * temperature field. The function fills a vector that contains average
       * field values over slices of the domain of same depth.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_temperature_averages(std::vector<double> &values) const;

      /**
       * Fill the argument with a set of lateral averages of the current
       * compositional fields.
       *
       * @param composition_index The index of the compositional field whose
       * matrix we want to assemble (0 <= composition_index < number of
       * compositional fields in this problem).
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_composition_averages(const unsigned int composition_index,
                               std::vector<double> &values) const;

      /**
       * Compute a lateral average of the current viscosity.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_viscosity_averages(std::vector<double> &values) const;

      /**
       * Compute a lateral average of the current velocity magnitude.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_velocity_magnitude_averages(std::vector<double> &values) const;

      /**
       * Compute a lateral average of the current sinking velocity.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_sinking_velocity_averages(std::vector<double> &values) const;

      /**
       * Compute a lateral average of the seismic shear wave speed: Vs.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_Vs_averages(std::vector<double> &values) const;

      /**
       * Compute a lateral average of the seismic pressure wave speed: Vp.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_Vp_averages(std::vector<double> &values) const;

      /**
       * Compute a lateral average of the vertical heat flux, with the sign
       * convention of positive heat flux when it flows upwards.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_vertical_heat_flux_averages(std::vector<double> &values) const;

      /**
       * Compute a lateral average of the vertical mass flux. Note that
       * while get_vertical_heat_flux_averages() computes the average
       * vertical heat flux (positive or negative), this function computes
       * the average of the total mass flux through a certain layer
       * (both down- and upward motion counted positively).
       *
       * @param values The output vector of laterally averaged values.
       */
      void
      get_vertical_mass_flux_averages(std::vector<double> &values) const;
  };
}


#endif
