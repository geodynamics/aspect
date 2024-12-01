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


#include <aspect/postprocess/visualization/strain_rate.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      StrainRate<dim>::
      StrainRate ()
        :
        DataPostprocessorScalar<dim> ("strain_rate",
                                      update_gradients | update_quadrature_points),
        Interface<dim>("1/s")
      {}



      template <int dim>
      void
      StrainRate<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,
                ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,
                ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = input_data.solution_gradients[q][d];

            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
            computed_quantities[q](0) = std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.));
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StrainRate,
                                                  "strain rate",
                                                  "A visualization output object that generates output "
                                                  "for the norm of the strain rate, i.e., for the quantity "
                                                  "$\\sqrt{\\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)}$ "
                                                  "in the incompressible case and "
                                                  "$\\sqrt{[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]}$ "
                                                  "in the compressible case."
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
                                                  "Physical units: \\si{\\per\\second}.")
    }
  }
}
