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

#ifndef _aspect_material_model_equation_of_state_linearized_incompressible_h
#define _aspect_material_model_equation_of_state_linearized_incompressible_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      using namespace dealii;

      /**
       * A simplified, incompressible equation of state where the density depends linearly
       * on temperature and composition, using the equation
       * $\rho(p,T,\mathfrak c) = \left(1-\alpha (T-T_0)\right)\rho_0 + \sum_i \Delta\rho_i \; \mathfrak c_i.$ "
       * There is no pressure-dependence of the density, and all other material properties
       * relating to the equation of state are assumed to be constant and identical for each
       * composition.
       */
      template <int dim>
      class LinearizedIncompressible
      {
        public:
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
           */
          bool is_compressible () const;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          static
          void
          set_number_of_compositions (const unsigned int n_comp);

          static unsigned int number_of_compositions;

        private:
          double reference_rho;
          double reference_T;
          double thermal_alpha;
          double reference_specific_heat;
          std::vector<double> compositional_delta_rhos;
      };
    }
  }
}

#endif
