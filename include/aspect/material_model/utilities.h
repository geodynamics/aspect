/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_utilities_h
#define _aspect_material_model_utilities_h

#include <aspect/global.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * An namespace in which we define utility functions that
     * might be used in many different places in the material
     * model to prevent code duplication.
     */
    namespace MaterialUtilities
    {
      /**
       * For multicomponent material models: Given a vector of of compositional
       * fields of length N, this function returns a vector of volume fractions
       * of length N+1, corresponding to the volume fraction of a ``background
       * material'' as the first entry, and volume fractions for each of the input
       * fields as the following entries. The returned vector will sum to one.
       * If the sum of the compositional_fields is greater than
       * one, we assume that there is no background mantle (i.e., that field value
       * is zero). Otherwise, the difference between the sum of the compositional
       * fields and 1.0 is assumed to be the amount of background mantle.
       * Optionally, one can input a component mask that determines which of the
       * compositional fields to use during the computation (e.g. because
       * some fields contain non-volumetric quantities like strain,
       * porosity, or trace elements). By default, all fields are included.
       */
      std::vector<double>
      compute_volume_fractions(const std::vector<double> &compositional_fields,
                               const ComponentMask &field_mask = ComponentMask());



      /**
       * For multicomponent material models:
       * Enumeration for selecting which averaging scheme to use when
       * averaging the properties of different compositional fields.
       * Select between harmonic, arithmetic, geometric, and
       * maximum_composition. The max composition scheme simply uses the
       * viscosity of whichever field has the highest volume fraction.
       */
      enum CompositionalAveragingOperation
      {
        harmonic,
        arithmetic,
        geometric,
        maximum_composition
      };



      /**
       * Read the compositional averaging operation from the parameter file,
       * using the parameter name given in @p parameter_name, and return the
       * enum that corresponds to this operation.
       */
      CompositionalAveragingOperation
      parse_compositional_averaging_operation (const std::string &parameter_name,
                                               const ParameterHandler &prm);



      /**
       * For multicomponent material models:
       * Material models compute output quantities such as the viscosity, the
       * density, etc. For some models, these values depend strongly on the
       * composition, and more than one compositional field might have nonzero
       * values at a given quadrature point. This means that properties have to
       * be averaged based on the fractions of each compositional field present.
       * This function performs this type of averaging. The averaging is based
       * on the choice in @p average_type. Averaging is conducted over the
       * compositional fields given in @p volume_fractions. This means that
       * @p volume_fractions and @p parameter_values have to have the same size,
       * which would typically be the number of compositional fields used in the
       * simulation (with the potential addition of a background field, in case
       * the composition does not add up to 1). However, one might not want to
       * average over all fields, as in some cases compositional fields do not
       * represent a rock type, but other tracked quantities like the finite
       * strain, so the implementation is independent of the number of entries in
       * @p volume_fractions.
       */
      double average_value (const std::vector<double> &volume_fractions,
                            const std::vector<double> &parameter_values,
                            const CompositionalAveragingOperation &average_type);



      /**
       * A data structure with all inputs for the
       * MaterialModel::Interface::compute_drucker_prager_yielding() method.
       */
      struct DruckerPragerInputs
      {
        /**
        * Constructor. Initializes the various variables of this
        * structure with the input values.
        */
        DruckerPragerInputs(const double cohesion,
                            const double friction_angle,
                            const double pressure,
                            const double effective_strain_rate);

        const double cohesion;
        const double friction_angle;
        const double pressure;
        const double effective_strain_rate;
      };



      /**
       * A data structure with all outputs computed by the
       * MaterialModel::Interface::compute_drucker_prager_yielding() method.
       */
      struct DruckerPragerOutputs
      {
        /**
        * Constructor. Initializes the various variables of this
        * structure to NaNs.
        */
        DruckerPragerOutputs();

        double yield_strength;
        double plastic_viscosity;
        double viscosity_pressure_derivative;
      };



      /**
       * For material models with plasticity:
       * Function to compute the material properties in @p out given the
       * inputs in @p in according to the the Drucker-Prager yield criterion.
       */
      template <int dim>
      void
      compute_drucker_prager_yielding (const DruckerPragerInputs &in,
                                       DruckerPragerOutputs &out);

    }
  }
}


#endif
