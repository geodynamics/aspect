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


#include <aspect/postprocess/visualization/friction_heating.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      FrictionHeating<dim>::
      FrictionHeating ()
        :
        DataPostprocessorScalar<dim> ("friction_heating",
                                      update_values | update_gradients | update_q_points)
      {}



      template <int dim>
      void
      FrictionHeating<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (duh[0].size() == this->introspection().n_components,          ExcInternalError());

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_quadrature_points,
                                                                       this->n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_quadrature_points,
                                                                         this->n_compositional_fields());

        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];

          }

        this->get_material_model().evaluate(in, out);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const SymmetricTensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>()
                 :
                 in.strain_rate[q]);

            computed_quantities[q](0) = 2 * out.viscosities[q] *
                                        compressible_strain_rate * compressible_strain_rate;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(FrictionHeating,
                                                  "friction heating",
                                                  "A visualization output object that generates output "
                                                  "for the amount of friction heating often referred "
                                                  "to as $\\tau:\\epsilon$. More concisely, in the "
                                                  "incompressible case, the quantity that is output "
                                                  "is defined as "
                                                  "$\\eta \\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)$ "
                                                  "where $\\eta$ is itself a function of temperature, "
                                                  "pressure and strain rate. In the compressible case, "
                                                  "the quantity that's computed is "
                                                  "$\\eta [\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]$.")
    }
  }
}
