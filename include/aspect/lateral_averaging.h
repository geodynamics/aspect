/*
  Copyright (C) 2011-2015 by the authors of the ASPECT code.

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


#ifndef __aspect__lateral_averaging_h
#define __aspect__lateral_averaging_h

#include <aspect/simulator_access.h>

namespace aspect
{
  using namespace dealii;

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
       * @param composition_index The index of the compositional field
       * to laterally average.
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
       * Compute a lateral average of the heat flux, with the sign
       * convention of positive heat flux when it flows upwards.
       *
       * @param values The output vector of laterally averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_vertical_heat_flux_averages(std::vector<double> &values) const;

    private:

      /**
       * Internal routine to compute the depth average of a certain quantitiy.
       *
       * The functor @p fctr must be an object of a user defined type that can
       * be arbitrary but has to satisfy certain requirements. In essence,
       * this class type needs to implement the following interface of member
       * functions:
       * @code
       * template <int dim>
       * class Functor
       * {
       *   public:
       *     // operator() will have @p in and @p out filled out if @p true
       *     bool need_material_properties() const;
       *
       *     // called once at the beginning with the number of quadrature points
       *     void setup(const unsigned int q_points);
       *
       *     // Fill @p output for each quadrature point.
       *     // This takes material model inputs and outputs (which are filled
       *     // if need_material_properties() == true), an initialized FEValues
       *     // object for a cell, and the current solution vector as inputs.
       *     // It then evaluates the desired quantity and puts the results in
       *     // the output vector, which is q_points long.
       *     void operator()(const MaterialModel::MaterialModelInputs<dim> &in,
       *                     const MaterialModel::MaterialModelOutputs<dim> &out,
       *                     FEValues<dim> &fe_values,
       *                     const LinearAlgebra::BlockVector &solution,
       *                     std::vector<double> &output);
       * };
       * @endcode
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       * @param fctr Instance of a class satisfying the signature above.
       */
      template<class FUNCTOR>
      void compute_lateral_average(std::vector<double> &values,
                                   FUNCTOR &fctr) const;

      /**
       * Compute a depth average of the current temperature/composition. The
       * function fills a vector that contains average
       * temperatures/compositions over slices of the domain of same depth.
       *
       * @param field  Extractor for temperature or compositional field to average.
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void get_field_averages(const FEValuesExtractors::Scalar &field,
                              std::vector<double> &values) const;

  };
}


#endif
