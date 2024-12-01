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
#include <aspect/boundary_fluid_pressure/interface.h>

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
        const double permeability = 1.0 + Utilities::fixed_power<2>(0.5);
        return permeability / 2.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {

        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            out.viscosities[i] = 3.0/4.0;
            out.densities[i] = 4.0;
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
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                melt_out->compaction_viscosities[i] = 1.0;
                melt_out->fluid_viscosities[i]= 2.0;
                melt_out->permeabilities[i]= 1.0 + Utilities::fixed_power<2>(in.position[i][1]);
                melt_out->fluid_densities[i]= 2.0;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
              }
          }
      }

  };

  template <int dim>
  class PressureBdry:

    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id boundary_indicator,
        const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim>> &normal_vectors,
        std::vector<double> &output
      ) const
      {
        for (unsigned int q=0; q<output.size(); ++q)
          {
            Tensor<1,dim> gradient;
            gradient[0] = 0.0;
            gradient[1] = -1.0;
            output[q] = gradient * normal_vectors[q];
          }
      }
  };

}



// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(MeltMaterial,
                                 "MeltMaterial3",
                                 "")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")
}
