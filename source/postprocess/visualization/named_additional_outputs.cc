/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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
        DataPostprocessor<dim> (),
        Interface<dim>("")  // physical units depend on run-time parameters
      {}



      template <int dim>
      void
      NamedAdditionalOutputs<dim>::
      initialize ()
      {
        MaterialModel::MaterialModelOutputs<dim> out(0,
                                                     this->n_compositional_fields());
        this->get_material_model().create_additional_named_outputs(out);

        AssertThrow(out.additional_outputs.size() > 0,
                    ExcMessage("You activated the postprocessor <named additional outputs>, but there are no additional outputs "
                               "provided by the material model. Either remove the postprocessor, or check why no output is provided."));

        for (unsigned int k=0; k<out.additional_outputs.size(); ++k)
          {
            const MaterialModel::NamedAdditionalMaterialOutputs<dim> *result
              = dynamic_cast<const MaterialModel::NamedAdditionalMaterialOutputs<dim> *> (out.additional_outputs[k].get());

            if (result)
              {
                std::vector<std::string> names = result->get_names();

                for (const auto &name : names)
                  property_names.push_back(name);
              }
          }

        AssertThrow(property_names.size() > 0,
                    ExcMessage("You activated the postprocessor <named additional outputs>, but none of the additional outputs "
                               "provided by the material model are named outputs. Either remove the postprocessor, or check why no "
                               "named output is provided."));

        for (auto &property_name : property_names)
          std::replace(property_name.begin(),property_name.end(),' ', '_');
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
        return update_gradients | update_values  | update_quadrature_points;
      }



      template <int dim>
      void
      NamedAdditionalOutputs<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,
                ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        this->get_material_model().create_additional_named_outputs(out);
        this->get_material_model().evaluate(in, out);

        unsigned int field_index = 0;
        for (unsigned int k=0; k<out.additional_outputs.size(); ++k)
          {
            const MaterialModel::NamedAdditionalMaterialOutputs<dim> *result
              = dynamic_cast<const MaterialModel::NamedAdditionalMaterialOutputs<dim> *> (out.additional_outputs[k].get());

            if (result)
              {
                std::vector<double> outputs(n_quadrature_points);
                for (unsigned int i=0; i<result->get_names().size(); ++i, ++field_index)
                  {
                    outputs = result->get_nth_output(i);

                    for (unsigned int q=0; q<n_quadrature_points; ++q)
                      computed_quantities[q][field_index] = outputs[q];
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
                                                  "in fact does not contain any named outputs at all."
                                                  "\n\n"
                                                  "Physical units: Various, depending on what is being output.")
    }
  }
}
