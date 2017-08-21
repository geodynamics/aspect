/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#include <aspect/heating_model/shear_heating_with_melt.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>
#include <deal.II/numerics/fe_field_function.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ShearHeatingMelt<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.strain_rate.size(),
             ExcMessage ("The shear heating plugin needs the strain rate!"));

      Assert(this->introspection().compositional_name_exists("porosity"),
             ExcMessage("Heating model shear heating with melt only works if there "
                        "is a compositional field called porosity."));

      // get the melt velocity from the solution vector
      std::vector<Tensor<1,dim> > melt_velocity (material_model_inputs.position.size());

      if (material_model_inputs.current_cell.state() == IteratorState::valid
          && this->get_timestep_number() > 0)
        {
          // we have to create a long vector, because that is the only way to extract the velocities
          // from the solution vector
          std::vector<std::vector<double> > melt_velocity_vector (dim, std::vector<double>(material_model_inputs.position.size()));
          // Prepare the field function
          Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector>
          fe_value(this->get_dof_handler(), this->get_solution(), this->get_mapping());

          fe_value.set_active_cell(material_model_inputs.current_cell);

          for (unsigned int d=0; d<dim; ++d)
            fe_value.value_list(material_model_inputs.position,
                                melt_velocity_vector[d],
                                this->introspection().variable("fluid velocity").first_component_index+d);

          for (unsigned int i=0; i<material_model_inputs.position.size(); ++i)
            for (unsigned int d=0; d<dim; ++d)
              melt_velocity[i][d] = melt_velocity_vector[d][i];
        }

      const MaterialModel::MeltOutputs<dim> *melt_outputs = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
      Assert(melt_outputs != NULL, ExcMessage("Need MeltOutputs from the material model for shear heating with melt."));

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          const double porosity = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("porosity")];

          if (porosity >= this->get_melt_handler().melt_transport_threshold)
            heating_model_outputs.heating_source_terms[q] = melt_outputs->compaction_viscosities[q]
                                                            * pow(trace(material_model_inputs.strain_rate[q]),2)
                                                            +
                                                            (melt_outputs->permeabilities[q] > 0
                                                             ?
                                                             melt_outputs->fluid_viscosities[q] * porosity * porosity
                                                             / melt_outputs->permeabilities[q]
                                                             * (melt_velocity[q] - material_model_inputs.velocity[q])
                                                             * (melt_velocity[q] - material_model_inputs.velocity[q])
                                                             :
                                                             0.0);
          else
            heating_model_outputs.heating_source_terms[q] = 0.0;

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }


    template <int dim>
    void
    ShearHeatingMelt<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output) const
    {
      MeltHandler<dim>::create_material_model_outputs(output);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeatingMelt,
                                  "shear heating with melt",
                                  "Implementation of a standard model for shear heating "
                                  "of migrating melt, including bulk (compression) heating "
                                  "$\\xi \\left( \\nabla \\cdot \\mathbf u_s \\right)^2 $ "
                                  "and heating due to melt segregation "
                                  "$\\frac{\\eta_f \\phi^2}{k} \\left( \\mathbf u_f - \\mathbf u_s \\right)^2 $. "
                                  "For full shear heating, "
                                  "this has to be used in combination with the heating model "
                                  "shear heating to also include shear heating for "
                                  "the solid part.")
  }
}
