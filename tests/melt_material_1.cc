/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>


namespace aspect
{
  template <int dim>
  class MeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double porosity = 0.01;
        const double permeability = 1.0 * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);
        return permeability / 0.1;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {

        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            out.viscosities[i] = 1.0;
            out.densities[i] = 1.0 + (in.composition[i][0]);
            out.thermal_expansion_coefficients[i] = 1.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim>>();

        if (melt_out != nullptr)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                const double porosity = in.composition[i][porosity_idx];

                melt_out->compaction_viscosities[i] = 100.0;
                melt_out->fluid_viscosities[i]=0.1;
                melt_out->permeabilities[i]=1.0 * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);
                melt_out->fluid_densities[i]=.1;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
              }
          }
      }

  };


}



// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(MeltMaterial,
                                 "MeltMaterial",
                                 "")
}
