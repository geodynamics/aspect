/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
      void
      SurfaceStress<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == SymmetricTensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

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
            const SymmetricTensor<2,dim> stress = -2.*eta*deviatoric_strain_rate +
                                                  in.pressure[q] * unit_symmetric_tensor<dim>();
            for (unsigned int i=0; i<SymmetricTensor<2,dim>::n_independent_components; ++i)
              computed_quantities[q](i) = stress[stress.unrolled_to_component_indices(i)];
          }

        // average the values if requested
        const auto &viz = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Visualization<dim> >();
        if (!viz.output_pointwise_stress_and_strain())
          average_quantities(computed_quantities);

      }


      template <int dim>
      std::vector<std::string>
      SurfaceStress<dim>::get_names () const
      {
        std::vector<std::string> names;
        switch (dim)
          {
            case 2:
              names.emplace_back("stress_xx");
              names.emplace_back("stress_yy");
              names.emplace_back("stress_xy");
              break;

            case 3:
              names.emplace_back("stress_xx");
              names.emplace_back("stress_yy");
              names.emplace_back("stress_zz");
              names.emplace_back("stress_xy");
              names.emplace_back("stress_xz");
              names.emplace_back("stress_yz");
              break;

            default:
              Assert (false, ExcNotImplemented());
          }

        return names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      SurfaceStress<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (SymmetricTensor<2,dim>::n_independent_components,
           DataComponentInterpretation::component_is_scalar);
      }



      template <int dim>
      UpdateFlags
      SurfaceStress<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_values | update_quadrature_points;
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
                                                  "in the compressible case. Note that the convention of positive "
                                                  "compressive stress is followed. ")
    }
  }
}
