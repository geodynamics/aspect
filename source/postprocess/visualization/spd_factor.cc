/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/spd_factor.h>

#include <aspect/newton.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SPDFactor<dim>::
	  SPDFactor ()
        :
        DataPostprocessorScalar<dim> ("spd_factor",
                                      update_values | update_gradients | update_q_points)
      {}



      template <int dim>
      void
	  SPDFactor<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &solution_values,
                                         const std::vector<std::vector<Tensor<1,dim> > > &solution_gradients,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (solution_values[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (solution_gradients[0].size() == this->introspection().n_components,          ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = solution_gradients[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=solution_values[q][this->introspection().component_indices.pressure];
            in.temperature[q]=solution_values[q][this->introspection().component_indices.temperature];
            for (unsigned int d = 0; d < dim; ++d)
              {
                in.velocity[q][d]=solution_values[q][this->introspection().component_indices.velocities[d]];
                in.pressure_gradient[q][d] = solution_gradients[q][this->introspection().component_indices.pressure][d];
              }

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = solution_values[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().evaluate(in, out);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
        {
        	const double eta = out.viscosities[q];
        	const MaterialModel::MaterialModelDerivatives<dim> *derivatives = out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();;
        	AssertThrow(derivatives != NULL, ExcMessage ("Error: The newton method requires the derivatives"));

        	const SymmetricTensor<2,dim> viscosity_derivative_wrt_strain_rate = derivatives->viscosity_derivative_wrt_strain_rate[q];
        	const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];

          computed_quantities[q](0) = Utilities::compute_spd_factor<dim>(strain_rate, viscosity_derivative_wrt_strain_rate, eta);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SPDFactor,
                                                  "spd factor",
                                                  "A visualization output object that generates output "
                                                  "for the viscosity.")
    }
  }
}
