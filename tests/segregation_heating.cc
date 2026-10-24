// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/gravity_model/interface.h>
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
        const std::shared_ptr<const MaterialModel::MeltOutputs<dim>>
        melt_outputs = material_model_outputs.template get_additional_output_object<MaterialModel::MeltOutputs<dim>>();
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
