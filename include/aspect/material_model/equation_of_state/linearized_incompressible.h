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
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
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
          /**
           * A function that computes the output of the equation of state @p out
           * for all compositions, given the inputs in @p in and an index @p q.
           * More specifically, the inputs structure MaterialModelInputs contains
           * a number of vectors, one for each input property, containing the
           * material model inputs at a number of different locations. The equation
           * of state model only evaluates one location (for each function call),
           * and the index q determines which entry of each of these vectors of
           * inputs is used.
           * Using these inputs, the equation of state outputs are computed
           * individually for each composition given in the inputs, and the outputs
           * structure is filled with the equation of state outputs, which contain a
           * number of vectors, one for each output property, with each vector
           * containing the separate outputs for each composition.
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
          declare_parameters (ParameterHandler &prm,
                              const unsigned int n_compositions = 0);

          /**
           * Read the parameters this class declares from the parameter file.
           * The optional parameter @p n_compositions determines the maximum
           * number of compositions the equation of state is set up with,
           * and should have the same value as the parameter with the same
           * name in the declare_parameters() function.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const unsigned int n_compositions = 0);


        private:
          /**
           * The reference density $\rho_0$ used in the computation of the density.
           */
          double reference_rho;

          /**
           * The reference temperature $T_0$ used in the computation of the density.
           */
          double reference_T;

          /**
           * The constant thermal expansivity $\alpha$ used in the computation of the density.
           */
          double thermal_alpha;

          /**
           * The constant specific heat.
           */
          double reference_specific_heat;

          /**
           * The maximum number of compositions that densities will be computed for, and,
           * accordingly, the maximum number of compositions the equation of state can be
           * evaluated with.
           */
          unsigned int maximum_number_of_compositions;

          /**
           * The density difference $\Delta\rho_i$ between each composition $\mathfrak c_i$
           * and the background.
           */
          std::vector<double> compositional_delta_rhos;
      };
    }
  }
}

#endif
