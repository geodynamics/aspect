/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/surface_stress.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/material_model/viscoelastic.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SurfaceStress<dim>::
      SurfaceStress ()
        :
        DataPostprocessorTensor<dim> ("surface_stress",
                                      update_values | update_gradients | update_quadrature_points),
        Interface<dim>("Pa")
      {}



      template <int dim>
      void
      SurfaceStress<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == Tensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        this->get_material_model().create_additional_named_outputs(out);

        // We do not need to compute anything but the viscosity
        in.requested_properties = MaterialModel::MaterialProperties::viscosity | MaterialModel::MaterialProperties::additional_outputs;

        // Compute the viscosity...
        this->get_material_model().evaluate(in, out);

        // ...and use it to compute the stresses
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
            const SymmetricTensor<2,dim> deviatoric_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            const double eta = out.viscosities[q];

            // Compressive stress is positive in geoscience applications
            SymmetricTensor<2,dim> stress = -2.*eta*deviatoric_strain_rate +
                                            in.pressure[q] * unit_symmetric_tensor<dim>();

            // Add elastic stresses if existent
            if (this->get_parameters().enable_elasticity == true)
              {
                SymmetricTensor<2, dim> stress_0, stress_old;

                stress_0[0][0] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx")];
                stress_0[1][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy")];
                stress_0[0][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy")];

                stress_old[0][0] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx_old")];
                stress_old[1][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy_old")];
                stress_old[0][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy_old")];

                if (dim == 3)
                  {
                    stress_0[2][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz")];
                    stress_0[0][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz")];
                    stress_0[1][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz")];

                    stress_old[2][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz_old")];
                    stress_old[0][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz_old")];
                    stress_old[1][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz_old")];
                  }

                const MaterialModel::ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<MaterialModel::ElasticAdditionalOutputs<dim>>();

                const double shear_modulus = elastic_out->elastic_shear_moduli[q];

                // Retrieve the timestep ratio dtc/dte and the elastic viscosity,
                // only two material models support elasticity.
                double timestep_ratio = 0.;
                double elastic_viscosity = 0.;
                if (Plugins::plugin_type_matches<MaterialModel::ViscoPlastic<dim>>(this->get_material_model()))
                  {
                    const MaterialModel::ViscoPlastic<dim> &vp = Plugins::get_plugin_as_type<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model());
                    elastic_viscosity = vp.get_elastic_viscosity(shear_modulus);
                    timestep_ratio = vp.get_timestep_ratio();
                  }
                else if (Plugins::plugin_type_matches<MaterialModel::Viscoelastic<dim>>(this->get_material_model()))
                  {
                    const MaterialModel::Viscoelastic<dim> &ve = Plugins::get_plugin_as_type<const MaterialModel::Viscoelastic<dim>>(this->get_material_model());
                    elastic_viscosity = ve.get_elastic_viscosity(shear_modulus);
                    timestep_ratio = ve.get_timestep_ratio();
                  }
                else
                  AssertThrow(false, ExcMessage("The stress component statistics postprocessor cannot be used with elasticity for material models other than ViscoPlastic and Viscoelastic."));

                // The total stress of timestep t.
                // Both eta and the elastic viscosity have been scaled with the timestep ratio.
                stress = in.pressure[q] * unit_symmetric_tensor<dim>() - (2. * eta * deviatoric_strain_rate + eta / elastic_viscosity * stress_0 + (1. - timestep_ratio) * (1. - eta / elastic_viscosity) * stress_old);
              }

            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int e=0; e<dim; ++e)
                computed_quantities[q][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
                  = stress[d][e];
          }

        // average the values if requested
        const auto &viz = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Visualization<dim>>();
        if (!viz.output_pointwise_stress_and_strain())
          average_quantities(computed_quantities);
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SurfaceStress,
                                                  "surface stress",
                                                  "A visualization output object that generates output "
                                                  "on the surface of the domain "
                                                  "for the 3 (in 2d) or 6 (in 3d) components of the stress "
                                                  "tensor, i.e., for the components of the tensor "
                                                  "$-2\\eta\\varepsilon(\\mathbf u)+pI$ "
                                                  "in the incompressible case and "
                                                  "$-2\\eta\\left[\\varepsilon(\\mathbf u)-"
                                                  "\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I\\right]+pI$ "
                                                  "in the compressible case. If elasticity is included, "
                                                  "its contribution is accounted for. Note that the convention of positive "
                                                  "compressive stress is followed."
                                                  "The stress outputted on the surface of the domain will equal "
                                                  "the stress on the surface of the volume output if the parameter "
                                                  "'Point-wise stress and strain' in the Visualization subsection "
                                                  "is set to true. "
                                                  "\n\n"
                                                  "Physical units: \\si{\\pascal}.")
    }
  }
}
