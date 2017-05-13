/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
      SPDfactor<dim>::
      SPDfactor ()
        :
        DataPostprocessorScalar<dim> ("spd_factor",
                                      update_values | update_gradients | update_q_points)
      {}



      template <int dim>
      void
      SPDfactor<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,          ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        in.position = input_data.evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = input_data.solution_gradients[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q] = input_data.solution_values[q][this->introspection().component_indices.pressure];
            in.temperature[q] = input_data.solution_values[q][this->introspection().component_indices.temperature];
            for (unsigned int d = 0; d < dim; ++d)
              {
                in.velocity[q][d] = input_data.solution_values[q][this->introspection().component_indices.velocities[d]];
                in.pressure_gradient[q][d] = input_data.solution_gradients[q][this->introspection().component_indices.pressure][d];
              }

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[c]];
          }

        std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> > mmd(new MaterialModel::MaterialModelDerivatives<dim> (n_quadrature_points));
        out.additional_outputs.push_back(mmd);
        this->get_material_model().evaluate(in, out);
        const MaterialModel::MaterialModelDerivatives<dim> *derivatives = out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();
        AssertThrow(derivatives != NULL, ExcMessage ("Error: The newton method requires the derivatives"));

        for (unsigned int q=0; q<n_quadrature_points; ++q)
        {
        	const double eta = out.viscosities[q];
        	const SymmetricTensor<2,dim> viscosity_derivative_wrt_strain_rate = derivatives->viscosity_derivative_wrt_strain_rate[q];
        	const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
        	computed_quantities[q](0) = Utilities::compute_spd_factor<dim>(eta, strain_rate, viscosity_derivative_wrt_strain_rate, 0.9);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SPDfactor,
                                                  "spd factor",
                                                  "A visualization output object that generates output "
                                                  "for the spd factor.")
    }
  }
}
