/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/surface_strain_rate_tensor.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      void
      SurfaceStrainRateTensor<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == SymmetricTensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = input_data.solution_gradients[q][d];

            const SymmetricTensor<2,dim> strain_rate = symmetrize(grad_u);
            const SymmetricTensor<2,dim> deviatoric_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            for (unsigned int i=0; i<SymmetricTensor<2,dim>::n_independent_components; ++i)
              computed_quantities[q](i) = deviatoric_strain_rate[deviatoric_strain_rate.unrolled_to_component_indices(i)];
          }

        // average the values if requested
        const auto &viz = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Visualization<dim> >();
        if (!viz.output_pointwise_stress_and_strain())
          average_quantities(computed_quantities);

      }



      template <int dim>
      std::vector<std::string>
      SurfaceStrainRateTensor<dim>::get_names () const
      {
        std::vector<std::string> names;
        switch (dim)
          {
            case 2:
              names.emplace_back("surface_strain_rate_tensor_xx");
              names.emplace_back("surface_strain_rate_tensor_yy");
              names.emplace_back("surface_strain_rate_tensor_xy");
              break;

            case 3:
              names.emplace_back("surface_strain_rate_tensor_xx");
              names.emplace_back("surface_strain_rate_tensor_yy");
              names.emplace_back("surface_strain_rate_tensor_zz");
              names.emplace_back("surface_strain_rate_tensor_xy");
              names.emplace_back("surface_strain_rate_tensor_xz");
              names.emplace_back("surface_strain_rate_tensor_yz");
              break;

            default:
              Assert (false, ExcNotImplemented());
          }

        return names;
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      SurfaceStrainRateTensor<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (SymmetricTensor<2,dim>::n_independent_components,
           DataComponentInterpretation::component_is_scalar);
      }



      template <int dim>
      UpdateFlags
      SurfaceStrainRateTensor<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_quadrature_points;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SurfaceStrainRateTensor,
                                                  "surface strain rate tensor",
                                                  "A visualization output object that generates output "
                                                  "on the surface of the domain "
                                                  "for the 3 (in 2d) or 6 (in 3d) components of the strain rate "
                                                  "tensor, i.e., for the components of the tensor "
                                                  "$\\varepsilon(\\mathbf u)$ "
                                                  "in the incompressible case and "
                                                  "$\\varepsilon(\\mathbf u)-"
                                                  "\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I$ "
                                                  "in the compressible case.")
    }
  }
}
