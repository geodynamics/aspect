/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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
#include <aspect/material_model/interface.h>



namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A data structure containing output fields that can be filled by the
     * evaluate() function of an EquationOfState model. It contains those
     * properties of the MaterialModelOutputs that are connected to the
     * equation of state (or the thermodynamic properties). Output values
     * are computed separately for each composition and phase.
     *
     * Accordingly, the vectors are the values for each compositional field
     * and individual phase at one specific location.
     */
    template <int dim>
    struct EquationOfStateOutputs
    {
      /**
       * Constructor. Initialize the various arrays of this structure with the
       * given number of compositions and phases.
       *
       * @param n_individual_compositions_and_phases The number of vector quantities
       * for which input will be provided, and outputs should be filled. Note that this
       * number does not have to be the number of compositions, it can be smaller (if
       * some compositional fields do not represent volumetric compositions, but tracked
       * quantities like strain) or larger (if there is a background field, or some
       * compositions have several high-pressure phases). It is the responsibility
       * of the material model and equation of state object to interpret the
       * entries consistently.
       */
      EquationOfStateOutputs (const unsigned int n_individual_compositions_and_phases);

      /**
       * Density values for each composition and phase.
       */
      std::vector<double> densities;

      /**
       * Thermal expansion coefficients for each composition and phase. It is defined
       * as $\alpha = - \frac{1}{\rho} \frac{\partial\rho}{\partial T}$
       */
      std::vector<double> thermal_expansion_coefficients;

      /**
       * Specific heat for each composition and phase.
       */
      std::vector<double> specific_heat_capacities;

      /**
       * Compressibility for each composition and phase. The compressibility is defined
       * as $\kappa = \frac{1}{\rho} \frac{\partial\rho}{\partial p}$.
       */
      std::vector<double> compressibilities;

      /**
       * The product of the change of entropy $\Delta S$ at a phase transition
       * and the derivative of the phase function $X=X(p,T,\mathfrak c,\mathbf
       * x)$ with regard to pressure for each composition and phase.
       */
      std::vector<double> entropy_derivative_pressure;

      /**
       * The product of (minus) the change of entropy $-\Delta S$ at a phase
       * transition and the derivative of the phase function
       * $X=X(p,T,\mathfrak c,\mathbf x)$ with regard to temperature for
       * each composition and phase.
       */
      std::vector<double> entropy_derivative_temperature;
    };



    /**
     * This function takes the output of an equation of state @p eos_outputs_all_phases,
     * which contains the data for all compositions and all of their phases at the
     * current conditions and uses a PhaseFunction object @p phase_function to compute
     * the effective value of equation of state properties for each individual composition.
     * Essentially it computes which phase is stable at the current conditions
     * (described in @p phase_in) and fills the equation of state output @p eos_outputs with
     * the properties of these stable phases.
     * @p eos_outputs now only contains as many entries as volumetric compositions (potentially
     * plus one for the background field, if the input @p eos_outputs_all_phases contained
     * one and the phase function contained one as well).
     */
    template <int dim>
    void
    phase_average_equation_of_state_outputs(const EquationOfStateOutputs<dim> &eos_outputs_all_phases,
                                            const std::vector<double> &phase_function_values,
                                            const std::vector<unsigned int> &n_phase_transitions_per_composition,
                                            EquationOfStateOutputs<dim> &eos_outputs);
  }
}

#endif
