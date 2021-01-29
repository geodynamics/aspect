/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_equation_of_state_multicomponent_incompressible_h
#define _aspect_material_model_equation_of_state_multicomponent_incompressible_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      using namespace dealii;

      /**
       * An incompressible equation of state that is intended for use with multiple compositional
       * fields and potentially phases. For each material property, the user supplies a comma
       * delimited list of length N+P+1, where N is the number of compositional fields used in
       * the computation, P is the total number of phase transitions.
       * The first entry corresponds to the "background" (which is also why there are N+P+1 entries).
       *
       * If a single value is given, then all the compositional fields and phases are given
       * that value. Other lengths of lists are not allowed. For a given
       * compositional field and phase the material parameters are treated as constant,
       * except density, which varies linearly with temperature according to the equation:
       *
       * $\rho(p,T,\mathfrak c) = \left(1-\alpha_i (T-T_{\text{Sref}})\right) \rho_0(\mathfrak c_i).$
       *
       * $T_{\text{Sref}}$ is the temperature on the reference adiabat
       * at the local pressure T(p, S_{\text{ref}}). If adiabatic heating is not included,
       * this temperature is constant and equal to the Adiabatic surface temperature.
       * In this case, the isothermal compressibility is equal to zero.
       * If adiabatic heating is included, $T_{\text{Sref}}$ changes as a function of
       * pressure such that there is no density change along the reference isentrope
       * (i.e., the isentropic compressibility is equal to zero).
       *
       * This equation of state is pleasingly simple but not thermodynamically self-consistent,
       * because the form of the density equation implies a non-constant thermal expansivity,
       * and the thermal expansivity divided by the pressure has a temperature dependence
       * which implies a non-constant isobaric heat capacity. These inconsistencies are
       * common to most equations of state used in geodynamic codes, and ubiquitous
       * among incompressible models.
       */
      template <int dim>
      class MulticomponentIncompressible :  public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * A function that computes the output of the equation of state @p out
           * for all compositions and phases, given the inputs in @p in and an
           * index input_index that determines which entry of the vector of inputs is used.
           */
          void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                        const unsigned int input_index,
                        MaterialModel::EquationOfStateOutputs<dim> &out) const;

          /**
           * Return whether the model is compressible or not. Incompressibility
           * does not necessarily imply that the density is constant; rather, it
           * may still depend on temperature or pressure. In the current
           * context, compressibility means whether we should solve the continuity
           * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
           * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
           */
          bool is_compressible () const;

          /**
           * Declare the parameters this class takes through input files.
           * The optional parameter @p default_thermal_expansion determines
           * the default value of the thermal expansivity used in the
           * equation of state.
           */
          static
          void
          declare_parameters (ParameterHandler &prm,
                              const double default_thermal_expansion = 3.5e-5);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and will be checked against the parsed
           * parameters.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition =
                              std::shared_ptr<std::vector<unsigned int>>());

        private:
          /**
           * Vector of reference densities $\rho_0$ with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> densities;

          /**
           * Vector of thermal expansivities with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> thermal_expansivities;

          /**
           * Vector of specific heat capacities with one entry per composition and phase plus one
           * for the background field.
           */
          std::vector<double> specific_heats;

          /**
           * A vector that stores how many separate phases there are per compositional
           * field. In other words if there is a background field without phase transitions
           * and one more composition that has 2 phase transitions this vector would store {1,3}.
           */
          std::shared_ptr<std::vector<unsigned int> > n_phases_per_composition;
      };
    }
  }
}

#endif
