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
      using namespace dealii;

      /**
       * A compressible equation of state that is intended for use with multiple compositional
       * fields. For each material property, the user supplies a comma delimited list of
       * length N+1, where N is the number of compositional fields used in the computation.
       * The first entry corresponds to the "background" (which is also why there are N+1 entries).
       *
       * If a single value is given, then all the compositional fields are given
       * that value. Other lengths of lists are not allowed. The material parameters
       * for each compositional field are calculated self-consistently,
       * assuming a constant pressure derivative of the bulk modulus at the reference temperature
       * (i.e. a Murnaghan equation of state), a constant ratio of the
       * thermal expansivity and isothermal compressibility, and a constant isochoric
       * specific heat. This leads to the following expressions for the material
       * properties of each material:
       *
       * $\rho(p,T) = \rho_0 f^{1/K_T'}$
       * $C_p(p,T) = C_v + (\alpha T \frac{\alpha}{\beta} f^{-1-(1/K_T')} / \rho_0)$
       * $\alpha(p, T) = \alpha_0/f$
       * $\beta(p, T) = \beta_0/f$
       * $f = (1 + (p - \frac{\alpha}{\beta}(T-T_0)) K_T' \beta)$.
       *
       */
      template <int dim>
      class MulticomponentCompressible :  public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * A function that computes the output of the equation of state @p out
           * for all compositions, given the inputs in @p in and an index q that
           * determines which entry of the vector of inputs is used.
           */
          void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                        const unsigned int q,
                        MaterialModel::EquationOfStateOutputs<dim> &out) const;

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
           * The optional parameter @p n_compositions determines the maximum
           * number of compositions the equation of state is set up with,
           * and should have the same value as the parameter with the same
           * name in the declare_parameters() function.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * Vector for reference_densities, read from parameter file .
           */
          std::vector<double> reference_densities;

        private:
          /**
           * Vector for reference temperatures, read from parameter file .
           */
          std::vector<double> reference_temperatures;

          /**
           * Vector for reference compressibilities, read from parameter file.
           */
          std::vector<double> reference_isothermal_compressibilities;

          /**
           * Vector for isothermal bulk modulus pressure derivatives, read from parameter file.
           */
          std::vector<double> isothermal_bulk_modulus_pressure_derivatives;

          /**
           * Vector for reference thermal expansivities, read from parameter file.
           */
          std::vector<double> reference_thermal_expansivities;

          /**
           * Vector for isochoric specific heats, read from parameter file.
           */
          std::vector<double> isochoric_specific_heats;

      };
    }
  }
}

#endif
