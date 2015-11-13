/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__averaging_h
#define __aspect__averaging_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
    * An enum to define what kind of averaging operations are implemented.
    * These are:
    *
    * - No averaging, i.e., leave the values as they were provided by the
    * material model.
    *
    * - Arithmetic averaging: Set the values of each output quantity at
    * every quadrature point to \f[ \bar x = \frac 1Q \sum_{q=1}^Q x_q \f]
    * where $x_q$ are the values at the $Q$ quadrature points.
    *
    * - Harmonic averaging: Set the values of each output quantity at every
    * quadrature point to \f[ \bar x = \left(\frac 1Q \sum_{q=1}^Q
    * \frac{1}{x_q}\right)^{-1} \f] where $x_q$ are the values at the $Q$
    * quadrature points.
    *
    * - Geometric averaging: Set the values of each output quantity at
    * every quadrature point to \f[ \bar x = \left(\prod_{q=1}^Q
    * x_q\right)^{1/Q} \f] where $x_q$ are the values at the $Q$ quadrature
    * points.
    *
    * - Pick largest: Set the values of each output quantity at every
    * quadrature point to \f[ \bar x = \max_{1\le q\le Q} x_q \f] where
    * $x_q$ are the values at the $Q$ quadrature points.
    *
    * - Log average: Set the values of each output quantity at every
    * quadrature point to \f[ \bar x = {10}^{\frac 1Q \sum_{q=1}^Q \log_{10} x_q} \f]
    * where $x_q$ are the values at the $Q$ quadrature points.
    *
    * - NWD Arithmetic averaging: Set the values of each output quantity at
    * every quadrature point to \f[ \bar x = \frac {\sum_{q=1}^Q W_q * x_q} {\sum_{q=1}^Q W_q} \f]
    * where $x_q$ are the values and $w$ the weights at the $Q$ quadrature points.
    *
    * - NWD Harmonic averaging: Set the values of each output quantity at every
    * quadrature point to \f[ \bar x = \frac{\sum_{q=1}^Q w_q }{ \sum_{q=1}^Q \frac{w_q}{x_q}} \f]
    *  where $x_q$ are the values and $w$ the weights
    * at the $Q$ quadrature points.
    *
    * - NWD Geometric averaging: Set the values of each output quantity at
    * every quadrature point to \f[ \bar x = \frac{ \sum_{q=1}^Q w_q x_q}
    * {\sum_{q=1}^Q w_q} \f] where $x_q$ are the values and $w$ the weights
    *  at the $Q$ quadrature points.
    */
    enum AveragingOperation
    {
      none,
      arithmetic_average,
      harmonic_average,
      geometric_average,
      pick_largest,
      log_average,
      nwd_arithmetic_average,
      nwd_harmonic_average,
      nwd_geometric_average
    };

    /**
     * A material model that applies an average of the quadrature points in a cell  to
     * a ''base model'' chosen from any of the other available material models.
     * @ingroup MaterialModels
     */

    template <int dim>
    class Averaging : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        virtual
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const;
        /**
         * Method to declare parameters related to depth-dependent model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to depth-dependent model
         */
        virtual void
        parse_parameters (ParameterHandler &prm);

        /**
         * Method that indicates whether material is compressible. Depth dependent model is compressible
         * if and only if base model is compressible.
         */
        virtual bool is_compressible () const;

        /**
         * Method to calculate reference viscosity for the depth-dependent model. The reference
         * viscosity is determined by evaluating the depth-dependent part of the viscosity at
         * the mean depth of the model.
         */
        virtual double reference_viscosity () const;

        /**
         * Method to calculate reference density for the depth-dependent model. Because the depth-
         * dependent model deos not modify density, the reference density is equivalent to the
         * base model's reference density.
         */
        virtual double reference_density () const;



      private:

        /**
        * Parse a string representing one of the options returned by
        * get_averaging_operation_names(), and return the corresponding
        * AveragingOperation value.
        */
        AveragingOperation
        parse_averaging_operation_name (const std::string &s);

        /**
        * Given the averaging @p operation, a description of where the
        * quadrature points are located on the given cell, and a mapping,
        * perform this operation on all elements of the @p values structure.
        */
        void
        average (const AveragingOperation averaging_operation,
                 const std::vector<Point<dim> >    &position,
                 std::vector<double>           &values_out) const;
        /**
        * The bell shape limit variable stores the maximum extend of the bell shape for the Normalized
        * Weighed Distance (NWD) averages.
        */
        double bell_shape_limit;
        /**
        * The averaging operation variable stores the chosen averaging operation.
        */
        AveragingOperation averaging_operation;
        /**
         * Pointer to the material model used as the base model
         */
        std_cxx11::shared_ptr<MaterialModel::Interface<dim> > base_model;
    };
  }
}

#endif
