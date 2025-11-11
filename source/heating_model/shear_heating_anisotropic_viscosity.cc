/*
 Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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

#include <aspect/heating_model/shear_heating_anisotropic_viscosity.h>

#include <aspect/material_model/interface.h>
#include <aspect/material_model/additional_outputs/anisotropic_viscosity.h>
#include <aspect/heating_model/shear_heating.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      // Some material models provide dislocation viscosities and boundary area work fractions
      // as additional material outputs. If they are attached, use them.
      const std::shared_ptr<const ShearHeatingOutputs<dim>> shear_heating_out =
        material_model_outputs.template get_additional_output_object<ShearHeatingOutputs<dim>>();

      const std::shared_ptr<const MaterialModel::AnisotropicViscosity<dim>> anisotropic_viscosity =
        material_model_outputs.template get_additional_output_object<MaterialModel::AnisotropicViscosity<dim>>();

      Assert(anisotropic_viscosity != nullptr,
             ExcMessage("This heating plugin should only be used with material models that provide "
                        "an anisotropic viscosity tensor, but none was provided."));

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // If there is an anisotropic viscosity, use it to compute the correct stress
          const SymmetricTensor<2,dim> &directed_strain_rate = anisotropic_viscosity->stress_strain_directors[q]
                                                               * material_model_inputs.strain_rate[q];

          const SymmetricTensor<2,dim> stress =
            2 * material_model_outputs.viscosities[q] *
            (this->get_material_model().is_compressible()
             ?
             directed_strain_rate - 1./3. * trace(directed_strain_rate) * unit_symmetric_tensor<dim>()
             :
             directed_strain_rate);

          const SymmetricTensor<2,dim> deviatoric_strain_rate =
            (this->get_material_model().is_compressible()
             ?
             material_model_inputs.strain_rate[q]
             - 1./3. * trace(material_model_inputs.strain_rate[q]) * unit_symmetric_tensor<dim>()
             :
             material_model_inputs.strain_rate[q]);

          heating_model_outputs.heating_source_terms[q] = stress * deviatoric_strain_rate;

          // If shear heating work fractions are provided, reduce the
          // overall heating by this amount (which is assumed to be converted into other forms of energy)
          if (shear_heating_out != nullptr)
            heating_model_outputs.heating_source_terms[q] *= shear_heating_out->shear_heating_work_fractions[q];

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      const unsigned int n_points = material_model_outputs.viscosities.size();

      if (material_model_outputs.template has_additional_output_object<MaterialModel::AnisotropicViscosity<dim>>() == false)
        {
          material_model_outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeatingAnisotropicViscosity,
                                  "anisotropic shear heating",
                                  "Implementation of a standard model for shear heating extended for an "
                                  "anisotropic viscosity tensor. If the material model provides a stress-"
                                  "strain director tensor, then the strain-rate is multiplied with this "
                                  "tensor to compute the stress that is used when computing the shear heating.")
  }
}
