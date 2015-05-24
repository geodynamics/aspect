/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/shear_stress.h>
#include <aspect/simulator_access.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      void
      ShearStress<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == SymmetricTensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (duh[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        // collect input information to compute the viscosity at every evaluation point
        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];
            for (unsigned int d = 0; d < dim; ++d)
              {
                in.velocity[q][d]=uh[q][this->introspection().component_indices.velocities[d]];
                in.pressure_gradient[q][d] = duh[q][this->introspection().component_indices.pressure][d];
              }

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        // then do compute the viscosity...
        this->get_material_model().evaluate(in, out);

        // ...and use it to compute the stresses
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
            const SymmetricTensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            const double eta = out.viscosities[q];

            const SymmetricTensor<2,dim> shear_stress = 2*eta*compressible_strain_rate;
            for (unsigned int i=0; i<SymmetricTensor<2,dim>::n_independent_components; ++i)
              computed_quantities[q](i) = shear_stress[shear_stress.unrolled_to_component_indices(i)];
          }
      }


      template <int dim>
      std::vector<std::string>
      ShearStress<dim>::get_names () const
      {
        std::vector<std::string> names;
        switch (dim)
          {
            case 2:
              names.push_back ("shear_stress_xx");
              names.push_back ("shear_stress_yy");
              names.push_back ("shear_stress_xy");
              break;

            case 3:
              names.push_back ("shear_stress_xx");
              names.push_back ("shear_stress_yy");
              names.push_back ("shear_stress_zz");
              names.push_back ("shear_stress_xy");
              names.push_back ("shear_stress_xz");
              names.push_back ("shear_stress_yz");
              break;

            default:
              Assert (false, ExcNotImplemented());
          }

        return names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      ShearStress<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (SymmetricTensor<2,dim>::n_independent_components,
           DataComponentInterpretation::component_is_scalar);
      }



      template <int dim>
      UpdateFlags
      ShearStress<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_values | update_q_points;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(ShearStress,
                                                  "shear stress",
                                                  "A visualization output object that generates output "
                                                  "for the 3 (in 2d) or 6 (in 3d) components of the shear stress "
                                                  "tensor, i.e., for the components of the tensor "
                                                  "$2\\eta\\varepsilon(\\mathbf u)$ "
                                                  "in the incompressible case and "
                                                  "$2\\eta\\left[\\varepsilon(\\mathbf u)-"
                                                  "\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I\\right]$ "
                                                  "in the compressible case. The shear "
                                                  "stress differs from the full stress tensor "
                                                  "by the absence of the pressure.")
    }
  }
}
