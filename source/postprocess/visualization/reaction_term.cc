/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/reaction_term.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      ReactionTerm<dim>::
      ReactionTerm ()
        :
        DataPostprocessor<dim> ()
      {}

      template <int dim>
      std::vector<std::string>
      ReactionTerm<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;
        for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
          solution_names.push_back (this->introspection().name_for_compositional_index(c) + "_change");

        return solution_names;
      }

      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      ReactionTerm<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
        for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
          interpretation.push_back (DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }

      template <int dim>
      UpdateFlags
      ReactionTerm<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_q_points;
      }

      template <int dim>
      void
      ReactionTerm<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                         const std::vector<Point<dim> >                  &normals,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == this->n_compositional_fields(), ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());

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
            for (unsigned int i = 0; i < dim; ++i)
              in.velocity[q][i]=uh[q][this->introspection().component_indices.velocities[i]];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().evaluate(in, out);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          for (unsigned int i=0; i<this->n_compositional_fields(); i++)
            computed_quantities[q][i] = out.reaction_terms[q][i];
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(ReactionTerm,
                                                  "reaction term",
                                                  "A visualization output object that generates output "
                                                  "for the reaction terms of all compositional fields.")
    }
  }
}
