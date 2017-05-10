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


#include <aspect/postprocess/visualization/named_additional_outputs.h>
#include <aspect/material_model/interface.h>
#include <aspect/utilities.h>

#include <algorithm>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      NamedAdditionalOutputs<dim>::
      NamedAdditionalOutputs ()
        :
        DataPostprocessor<dim> ()
      {}


      template <int dim>
      void
      NamedAdditionalOutputs<dim>::
      initialize ()
      {
        MaterialModel::MaterialModelOutputs<dim> out(0,
                                                     this->n_compositional_fields());
        this->get_material_model().create_additional_named_outputs(out);

        for (unsigned int k=0; k<out.additional_outputs.size(); ++k)
          {
            const MaterialModel::NamedAdditionalMaterialOutputs<dim> *result
              = dynamic_cast<const MaterialModel::NamedAdditionalMaterialOutputs<dim> *> (out.additional_outputs[k].get());

            if (result)
              {
                std::vector<std::string> names = result->get_names();

                for (unsigned int i=0; i<names.size(); ++i)
                  property_names.push_back(names[i]);
              }
          }

        for (unsigned int i=0; i<property_names.size(); ++i)
          std::replace(property_names[i].begin(),property_names[i].end(),' ', '_');
      }



      template <int dim>
      std::vector<std::string>
      NamedAdditionalOutputs<dim>::
      get_names () const
      {
        return property_names;
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      NamedAdditionalOutputs<dim>::
      get_data_component_interpretation () const
      {
        return std::vector<DataComponentInterpretation::DataComponentInterpretation> (get_names().size(),
                                                                                      DataComponentInterpretation::component_is_scalar);
      }



      template <int dim>
      UpdateFlags
      NamedAdditionalOutputs<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_q_points;
      }



      template <int dim>
      void
      NamedAdditionalOutputs<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,
                ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        in.position = input_data.evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              {
                grad_u[d] = input_data.solution_gradients[q][d];
                in.velocity[q][d] = input_data.solution_values[q][this->introspection().component_indices.velocities[d]];
                in.pressure_gradient[q][d] = input_data.solution_gradients[q][this->introspection().component_indices.pressure][d];
              }

            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q] = input_data.solution_values[q][this->introspection().component_indices.pressure];
            in.temperature[q] = input_data.solution_values[q][this->introspection().component_indices.temperature];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().create_additional_named_outputs(out);
        this->get_material_model().evaluate(in, out);

        for (unsigned int k=0; k<out.additional_outputs.size(); ++k)
          {
            const MaterialModel::NamedAdditionalMaterialOutputs<dim> *result
              = dynamic_cast<const MaterialModel::NamedAdditionalMaterialOutputs<dim> *> (out.additional_outputs[k].get());

            if (result)
              {
                std::vector<double> outputs(n_quadrature_points);
                for (unsigned int i=0; i<get_names().size(); ++i)
                  {
                    outputs = result->get_nth_output(i);

                    for (unsigned int q=0; q<n_quadrature_points; ++q)
                      computed_quantities[q][i] = outputs[q];
                  }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(NamedAdditionalOutputs,
                                                  "named additional outputs",
                                                  "Some material models can compute quantities other than those "
                                                  "that typically appear in the equations that \\aspect{} solves "
                                                  "(such as the viscosity, density, etc). Examples of quantities "
                                                  "material models may be able to compute are seismic velocities, "
                                                  "or other quantities that can be derived from the state variables "
                                                  "and the material coefficients such as the stress or stress "
                                                  "anisotropies. These quantities are generically referred to as "
                                                  "`named outputs' because they are given an explicit name different "
                                                  "from the usual outputs of material models.\n\n"
                                                  "This visualization postprocessor outputs whatever quantities the "
                                                  "material model can compute. What quantities these are is specific "
                                                  "to the material model in use for a simulation, and for many models "
                                                  "in fact does not contain any named outputs at all.")
    }
  }
}
