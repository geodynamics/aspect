/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/melt.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/melt_interface.h>

#include <deal.II/numerics/data_out.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
       MeltPressure<dim>::
       MeltPressure ()
         :
         DataPostprocessorScalar<dim> ("p_c",
                                       update_values | update_gradients | update_q_points)
       {}



       template <int dim>
       void
       MeltPressure<dim>::
       compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                          const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                          const std::vector<std::vector<Tensor<2,dim> > > &,
                                          const std::vector<Point<dim> > &,
                                          const std::vector<Point<dim> >                  &evaluation_points,
                                          std::vector<Vector<double> >                    &computed_quantities) const
       {
         // p_f = (p_c - (1-phi) p_s ) / (phi-1)
         // or p_c if phi=1
         // melt velocity = v_f =  v_s - K_D (nabla p_f - rho_f g) / phi  or = 0

         const unsigned int n_quadrature_points = uh.size();
         Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
         Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
         Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());
         Assert (duh[0].size() == this->introspection().n_components,          ExcInternalError());

         typename MaterialModel::MeltInterface<dim>::MaterialModelInputs in(n_quadrature_points, this->n_compositional_fields());
         typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs out(n_quadrature_points, this->n_compositional_fields());

         const unsigned int por_idx = this->introspection().compositional_index_for_name("porosity");

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

         const typename MaterialModel::MeltInterface<dim> * melt_mat = dynamic_cast<const MaterialModel::MeltInterface<dim>*> (&this->get_material_model());
         AssertThrow(melt_mat != NULL, ExcMessage("Need MeltMaterial if include_melt_transport is on."));
         melt_mat->evaluate_with_melt(in, out);

         for (unsigned int q=0; q<n_quadrature_points; ++q)
           {
             double phi = in.composition[q][por_idx];

             double p_s = in.pressure[q];
             double p_f = uh[q][this->introspection().component_indices.compaction_pressure];
             double p_c = (1.0-phi) * (p_s - p_f);

             computed_quantities[q](0) = p_c;
           }
       }

      template <int dim>
      MeltVelocity<dim>::
      MeltVelocity ()
        :
        DataPostprocessorVector<dim> ("melt_velocity",
                                      update_values | update_gradients | update_q_points)
      {}



      template <int dim>
      void
      MeltVelocity<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        // p_f = (p_c - (1-phi) p_s ) / (phi-1)
        // or p_c if phi=1
        // melt velocity = v_f =  v_s - K_D (nabla p_f - rho_f g) / phi  or = 0

        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == dim,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (duh[0].size() == this->introspection().n_components,          ExcInternalError());

        typename MaterialModel::MeltInterface<dim>::MaterialModelInputs in(n_quadrature_points, this->n_compositional_fields());
        typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs out(n_quadrature_points, this->n_compositional_fields());

        const unsigned int por_idx = this->introspection().compositional_index_for_name("porosity");

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

        const typename MaterialModel::MeltInterface<dim> * melt_mat = dynamic_cast<const MaterialModel::MeltInterface<dim>*> (&this->get_material_model());
        AssertThrow(melt_mat != NULL, ExcMessage("Need MeltMaterial if include_melt_transport is on."));
        melt_mat->evaluate_with_melt(in, out);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            double phi = in.composition[q][por_idx];

            double p_s = in.pressure[q];
            double p_f = uh[q][this->introspection().component_indices.compaction_pressure];

            if (phi < 1e-7)
              for (unsigned int d=0; d<dim; ++d)
                computed_quantities[q](d) = uh[q][d];
            else
              {
                double K_D = out.permeabilities[q] / out.fluid_viscosities[q];
                Tensor<1,dim> grad_p_f = duh[q][this->introspection().component_indices.compaction_pressure];
                const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(evaluation_points[q]);

                // v_f =  v_s - K_D (nabla p_f - rho_f g) / phi
                Tensor<1,dim> correction = - K_D * (grad_p_f - out.fluid_densities[q] * gravity) / phi;

                for (unsigned int d=0; d<dim; ++d)
                  computed_quantities[q](d) = uh[q][d] + correction[d];
              }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltVelocity,
                                                  "melt velocity",
                                                  "A visualization output object that generates output "
                                                  "for melt velocity.")
          ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltPressure,
                                                      "melt pressure",
                                                      "A visualization output object that generates output "
                                                      "for melt pressure.")
    }
  }
}
