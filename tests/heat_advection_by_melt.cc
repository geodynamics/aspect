/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/melt_global.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

namespace aspect
{

  template <int dim>
  class PressureBdry:
    public BoundaryFluidPressure::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id /*boundary_indicator*/,
        const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim>> &normal_vectors,
        std::vector<double> &fluid_pressure_gradient_outputs
      ) const
      {
        const MaterialModel::MeltOutputs<dim> *melt_outputs = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim>>();
        Assert(melt_outputs != nullptr, ExcMessage("Need MeltOutputs from the material model for shear heating with melt."));

        for (unsigned int q=0; q<fluid_pressure_gradient_outputs.size(); ++q)
          {
            const double density = material_model_outputs.densities[q];
            const double melt_density = melt_outputs->fluid_densities[q];
            const double porosity = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("porosity")];

            const double bulk_density = (porosity * melt_density + (1.0 - porosity) * density);
            fluid_pressure_gradient_outputs[q] = this->get_gravity_model().gravity_vector(material_model_inputs.position[q]) * bulk_density * normal_vectors[q];
          }
      }

  };

  namespace MaterialModel
  {
    template <int dim>
    class Advection : public MaterialModel::MeltGlobal<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
    };


    template <int dim>
    void
    Advection<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MeltGlobal<dim>::evaluate(in, out);

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          out.densities[i] = 2.0;
        }

      if (melt_out != nullptr)
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              melt_out->fluid_densities[i] = 1.0;
              melt_out->permeabilities[i] = 1.0;
            }
        }



    }

  }

}

// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Advection,
                                   "melt global advection",
                                   "A simple material model that is like the "
                                   "'melt simple' model, but has a constant reaction term.")
  }
}
