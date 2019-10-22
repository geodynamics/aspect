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

#ifndef _aspect_material_model_equation_of_state_interface_h
#define _aspect_material_model_equation_of_state_interface_h

#include <aspect/global.h>
#include <aspect/material_model/utilities.h>



namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A data structure containing output fields that can be filled by the
     * evaluate() function of an EquationOfState model. It contains those
     * properties of the MaterialModelOutputs that are connected to the
     * equation of state (or the thermodynamic properties). Output values
     * are computed separately for each composition.
     *
     * Accordingly, the vectors are the values for each compositional field
     * at one specific location.
     */
    template <int dim>
    struct EquationOfStateOutputs
    {
      /**
       * Constructor. Initialize the various arrays of this structure with the
       * given number of compositions.
       *
       * @param n_compositions The number of vector quantities (in the order in which
       * the Introspection class reports them) for which input will be
       * provided, and outputs should be filled.
       */
      EquationOfStateOutputs (const unsigned int n_compositions);

      /**
       * Density values for each composition.
       */
      std::vector<double> densities;

      /**
       * Thermal expansion coefficients for each composition. It is defined
       * as $\alpha = - \frac{1}{\rho} \frac{\partial\rho}{\partial T}$
       */
      std::vector<double> thermal_expansion_coefficients;

      /**
       * Specific heat for each composition.
       */
      std::vector<double> specific_heat_capacities;

      /**
       * Compressibility for each composition. The compressibility is defined
       * as $\kappa = \frac{1}{\rho} \frac{\partial\rho}{\partial p}$.
       */
      std::vector<double> compressibilities;

      /**
       * The product of the change of entropy $\Delta S$ at a phase transition
       * and the derivative of the phase function $X=X(p,T,\mathfrak c,\mathbf
       * x)$ with regard to pressure for each composition.
       */
      std::vector<double> entropy_derivative_pressure;

      /**
       * The product of (minus) the change of entropy $-\Delta S$ at a phase
       * transition and the derivative of the phase function
       * $X=X(p,T,\mathfrak c,\mathbf x)$ with regard to temperature for
       * each composition.
       */
      std::vector<double> entropy_derivative_temperature;
    };

    template <int dim>
    void
    compute_equation_of_state_phase_transitions(const EquationOfStateOutputs<dim> &eos_outputs_all_phases,
                                                const MaterialUtilities::PhaseFunction<dim> &phase_function,
                                                MaterialUtilities::PhaseFunctionInputs<dim> &phase_in,
                                                EquationOfStateOutputs<dim> &eos_outputs);
  }
}

#endif
