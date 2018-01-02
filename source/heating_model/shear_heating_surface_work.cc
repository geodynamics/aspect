/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/heating_model/shear_heating_surface_work.h>
#include <aspect/material_model/damage_rheology.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ShearHeatingSurfaceWork<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.strain_rate.size(),
             ExcMessage ("The shear heating plugin needs the strain rate!"));

      const MaterialModel::DislocationViscosityOutputs<dim> *disl_viscosities_out =
        material_model_outputs.template get_additional_output<MaterialModel::DislocationViscosityOutputs<dim> >();

      Assert(disl_viscosities_out != 0,
             ExcMessage ("The shear heating with surface work plugin needs the dislocation viscosities "
                         "as additional material model output. Currently this is only implemented in the 'damage rheology' "
                         "material model."));

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          const SymmetricTensor<2,dim> compressible_strain_rate
            = (this->get_material_model().is_compressible()
               ?
               material_model_inputs.strain_rate[q] - 1./3. * trace(material_model_inputs.strain_rate[q]) * unit_symmetric_tensor<dim>()
               :
               material_model_inputs.strain_rate[q]);

          // This is the usual heating term
          heating_model_outputs.heating_source_terms[q] = 2.0 * material_model_outputs.viscosities[q] *
                                                          compressible_strain_rate * compressible_strain_rate;

          // This is the part of the energy that is used to decrease the grain size
          // Only the deformation that is accomodated by dislocation creep contributes
          // to the grain size reduction.

          heating_model_outputs.heating_source_terms[q] -= 2.0 * material_model_outputs.viscosities[q] *
                                                           disl_viscosities_out->boundary_area_change_work_fraction[q] *
                                                           material_model_outputs.viscosities[q] / disl_viscosities_out->dislocation_viscosities[q] *
                                                           compressible_strain_rate * compressible_strain_rate;

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }

    template <int dim>
    void
    ShearHeatingSurfaceWork<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output) const
    {
      this->get_material_model().create_additional_named_outputs(output);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeatingSurfaceWork,
                                  "shear heating surface work",
                                  "Implementation of a model for shear heating."
                                  "This model removes part of the standard model "
                                  "description. This part is assumed to increase "
                                  "the surface energy of the grain assemblage and "
                                  "therefore does not contribute to the generated heat."
                                  "This plugin is limited to the damage rheology "
                                  "material model.")
  }
}
