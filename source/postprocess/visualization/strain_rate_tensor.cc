/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/strain_rate_tensor.h>


#if DEAL_II_VERSION_GTE(9,0,0)
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      StrainRateTensor<dim>::
      StrainRateTensor ()
        :
        DataPostprocessorTensor<dim> ("strain_rate_tensor",
                                      update_gradients | update_quadrature_points)
      {}

      template <int dim>
      void
      StrainRateTensor<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert ((computed_quantities[0].size() == Tensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = input_data.solution_gradients[q][d];

            const SymmetricTensor<2,dim> strain_rate = symmetrize(grad_u);
            const Tensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            for (unsigned int i=0; i<Tensor<2,dim>::n_independent_components; ++i)
              computed_quantities[q](i)
                = compressible_strain_rate[compressible_strain_rate.unrolled_to_component_indices(i)];
          }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StrainRateTensor,
                                                  "strain rate tensor",
                                                  "A visualization output object that generates output "
                                                  "for the 4 (in 2d) or 9 (in 3d) components of the strain rate "
                                                  "tensor, i.e., for the components of the tensor "
                                                  "$\\varepsilon(\\mathbf u)$ "
                                                  "in the incompressible case and "
                                                  "$\\varepsilon(\\mathbf u)-"
                                                  "\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I$ "
                                                  "in the compressible case.")
    }
  }
}
#endif
