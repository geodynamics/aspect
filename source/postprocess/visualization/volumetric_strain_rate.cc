/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/volumetric_strain_rate.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VolumetricStrainRate<dim>::
      VolumetricStrainRate ()
        :
        DataPostprocessorScalar<dim> ("volumetric_strain_rate",
                                      update_gradients | update_q_points)
      {}


      template <int dim>
      void
      VolumetricStrainRate<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,          ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // extract the velocity gradients
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = input_data.solution_gradients[q][d];

            const double div_u = trace(grad_u);
            computed_quantities[q](0) = div_u;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VolumetricStrainRate,
                                                  "volumetric strain rate",
                                                  "A visualization output object that generates output "
                                                  "for the volumetric strain rate, i.e., for the quantity "
                                                  "$\\nabla\\cdot\\mathbf u = \\textrm{div}\\; \\mathbf u = "
                                                  "\\textrm{trace}\\; \\varepsilon(\\mathbf u)$. "
                                                  "This should be zero (in some average sense) in incompressible "
                                                  "convection models, but can be non-zero in compressible models and "
                                                  "models with melt transport.")
    }
  }
}
