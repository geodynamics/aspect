/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_equation_of_state_multicomponent_compressible_h
#define _aspect_material_model_equation_of_state_multicomponent_compressible_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      /**
       * A compressible equation of state that is intended for use with multiple compositional
       * fields and potentially phases.
       *
       * For models without phase transitions, each material property is specified as a comma
       * delimited list of length N+1 for the background and compositional fields (N) For models
       * with phase transitions, the list needs to contain each field name, including the background,
       * for a total of N+1 names, and for each of these names, specify the value for each phase
       * Therefore, the total number of values given is N+P+1, with P = sum(P_c) the total number of
       * phase transitions, summed over all phases. The format is background: value1|value2|...|valueP_1+1,
       * field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. If only one value is given,
       * then all fields/phases use the same value. Other lengths of lists are not allowed.
       *
       * If no phase transitions are included, the material parameters for each compositional field
       * are calculated self-consistently, assuming a constant pressure derivative of the isothermal
       * bulk modulus ($K_T'$) at the reference temperature (i.e. a Murnaghan equation of state),
       * a constant ratio of the thermal expansivity ($\alpha$) and isothermal compressibility ($\beta_T$),
       * and a constant isochoric specific heat $C_v$. This leads to the following expressions for the material
       * properties of each material:
       *
       * $\rho(p,T) = \rho_0 f^{1/K_T'}$
       * $C_p(p,T) = C_v + (\alpha T \frac{\alpha}{\beta} f^{-1-(1/K_T')} / \rho_0)$
       * $\alpha(p, T) = \alpha_0/f$
       * $\beta_T(p, T) = \beta_0/f$
       * $f = (1 + (p - \frac{\alpha}{\beta}(T-T_0)) K_T' \beta)$.
       *
       * where $\rho$ is the density and $C_p$ is the isobaric heat capacity.
       * $f$ is a scaling factor for $\alpha$ and $\beta_T$.
       *
       * Significantly, if phase transitions are included the formulation no longer self-consistently calculates
       * the second derivative properties (heat capacity, thermal expansivity, and compressibility), as they are
       * they are affected by the P-T-X dependence of the phase function.
       */
      template <int dim>
      class MulticomponentCompressible :  public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * A function that computes the output of the equation of state @p eos_outputs
           * for all compositions and phases, given the inputs in @p in and an
           * index q that determines which entry of the vector of inputs is used.
           */
          void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                        const unsigned int q,
                        MaterialModel::EquationOfStateOutputs<dim> &eos_outputs) const;

          /**
           * Return whether the model is compressible or not. Incompressibility
           * does not necessarily imply that the density is constant; rather, it
           * may still depend on temperature or pressure. In the current
           * context, compressibility means whether we should solve the continuity
           * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
           * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
           * This model is compressible.
           */
          bool is_compressible () const;

          /**
           * Declare the parameters this class takes through input files.
           * The optional parameter @p n_compositions determines the maximum
           * number of compositions the equation of state is set up with,
           * in other words, how many compositional fields influence the
           * density.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers, this is considered the number of phases
           * for each compositional field and will be checked against the parsed
           * parameters.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition = nullptr);

          /**
           * Vector of reference densities $\rho_0$ with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> reference_densities;

        private:
          /**
           * Vector of reference temperatures $T_0$ with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> reference_temperatures;

          /**
           * Vector of reference compressibilities, with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> reference_isothermal_compressibilities;

          /**
           * Vector of isothermal bulk modulus pressure derivatives with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> isothermal_bulk_modulus_pressure_derivatives;

          /**
           * Vector of reference thermal expansivities, with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> reference_thermal_expansivities;

          /**
           * Vector of isochoric specific heats, with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> isochoric_specific_heats;

          /**
           * Whether to enable the use of phase transitions, which currently breaks thermodynamic consistency
           */
          bool enable_phase_transitions;
      };
    }
  }
}

#endif
