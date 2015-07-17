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


#include <aspect/postprocess/visualization/heating.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/utilities.h>

#include <algorithm>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Heating<dim>::
      Heating ()
        :
        DataPostprocessor<dim> ()
      {}

      template <int dim>
      std::vector<std::string>
      Heating<dim>::
      get_names () const
      {
        std::vector<std::string> names = this->get_heating_model_manager().get_active_heating_model_names();

        // make the names valid names for output variables via DataOut
        for (unsigned int i=0; i<names.size(); ++i)
          std::replace(names[i].begin(), names[i].end(), ' ', '_');

        return names;
      }

      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      Heating<dim>::
      get_data_component_interpretation () const
      {
        return std::vector<DataComponentInterpretation::DataComponentInterpretation>
               (this->get_heating_model_manager().get_active_heating_model_names().size(),
                DataComponentInterpretation::component_is_scalar);
      }

      template <int dim>
      UpdateFlags
      Heating<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_q_points;
      }

      template <int dim>
      void
      Heating<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        const std::list<std_cxx11::shared_ptr<HeatingModel::Interface<dim> > > &heating_model_objects = this->get_heating_model_manager().get_active_heating_models();

        // we do not want to write any output if there are no heating models
        // used in the computation
        if (heating_model_objects.size() == 0)
          return;

        Assert (computed_quantities.size() == n_quadrature_points,
                ExcMessage("The length of the vector of quantities that are computed in the "
                           "postprocessor (" + dealii::Utilities::int_to_string(computed_quantities.size()) + ") has to match the "
                           "number of quadrature points (" + dealii::Utilities::int_to_string(n_quadrature_points) + ")!"));
        Assert (computed_quantities[0].size() == heating_model_objects.size(), ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components, ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());

        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (n_quadrature_points));

        HeatingModel::HeatingModelOutputs heating_model_outputs(n_quadrature_points, this->n_compositional_fields());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.temperature[q] = uh[q][this->introspection().component_indices.temperature];
            in.pressure[q]    = uh[q][this->introspection().component_indices.pressure];

            for (unsigned int d = 0; d < dim; ++d)
              {
                in.velocity[q][d]=uh[q][this->introspection().component_indices.velocities[d]];
                in.pressure_gradient[q][d] = duh[q][this->introspection().component_indices.pressure][d];
              }

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        in.position = evaluation_points;

        this->get_material_model().evaluate(in, out);

        unsigned int index = 0;
        for (typename std::list<std_cxx11::shared_ptr<HeatingModel::Interface<dim> > >::const_iterator
             heating_model = heating_model_objects.begin();
             heating_model != heating_model_objects.end(); ++heating_model, ++index)
          {
            (*heating_model)->evaluate(in, out, heating_model_outputs);

            for (unsigned int q=0; q<n_quadrature_points; ++q)
              computed_quantities[q][index] = heating_model_outputs.heating_source_terms[q];
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Heating,
                                                  "heating",
                                                  "A visualization output object that generates output "
                                                  "for all the heating terms used in the energy equation.")
    }
  }
}
