/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

        // We do not need to compute anything but the viscosity
        in.requested_properties = MaterialModel::MaterialProperties::viscosity;

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
                stress[0][0] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx")];
                stress[1][1] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy")];
                stress[0][1] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy")];

                if (dim == 3)
                  {
                    stress[2][2] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz")];
                    stress[0][2] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz")];
                    stress[1][2] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz")];
                  }
              }

            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int e=0; e<dim; ++e)
                computed_quantities[q][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
                  = stress[d][e];
          }

        // average the values if requested
        const auto &viz = this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::Visualization<dim>>();
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
                                                  "This postprocessor outputs the quantity computed herein as "
                                                  "a tensor, i.e., programs such as VisIt or Pararview can "
                                                  "visualize it as tensors represented by ellipses, not just "
                                                  "as individual fields. That said, you can also visualize "
                                                  "individual tensor components, by noting that the "
                                                  "components that are written to the output file correspond to "
                                                  "the tensor components $t_{xx}, t_{xy}, t_{yx}, t_{yy}$ (in 2d) "
                                                  "or  $t_{xx}, t_{xy}, t_{xz}, t_{yx}, t_{yy}, t_{yz}, t_{zx}, t_{zy}, "
                                                  "t_{zz}$ (in 3d) of a tensor $t$ in a Cartesian coordinate system. "
                                                  "Even though the tensor we output is symmetric, the output contains "
                                                  "all components of the tensor because that is what the file format "
                                                  "requires."
                                                  "\n\n"
                                                  "Physical units: \\si{\\pascal}.")
    }
  }
}
