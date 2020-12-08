/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/heating.h>
#include <aspect/heating_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/grid_tools.h>

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
        return update_gradients | update_values  | update_quadrature_points | update_JxW_values;
      }

      template <int dim>
      void
      Heating<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        const auto &heating_model_objects = this->get_heating_model_manager().get_active_heating_models();

        // we do not want to write any output if there are no heating models
        // used in the computation
        if (heating_model_objects.size() == 0)
          return;

        Assert (computed_quantities.size() == n_quadrature_points,
                ExcMessage("The length of the vector of quantities that are computed in the "
                           "postprocessor (" + dealii::Utilities::int_to_string(computed_quantities.size()) + ") has to match the "
                           "number of quadrature points (" + dealii::Utilities::int_to_string(n_quadrature_points) + ")!"));
        Assert (computed_quantities[0].size() == heating_model_objects.size(), ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components, ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data, this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());

        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (n_quadrature_points));

        this->get_heating_model_manager().create_additional_material_model_inputs_and_outputs(in, out);
        HeatingModel::HeatingModelOutputs heating_model_outputs(n_quadrature_points, this->n_compositional_fields());

        // we need the cell as input for the material model because some heating models
        // want to access the solution vector.
#if DEAL_II_VERSION_GTE(9,3,0)
        in.current_cell = input_data.template get_cell<dim>();
#else
        in.current_cell = input_data.template get_cell<DoFHandler<dim> > ();
#endif

        // we need an fevalues object to get the melt velocities
        std::vector<Point<dim> > quadrature_points(n_quadrature_points);
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          quadrature_points[q] = this->get_mapping().transform_real_to_unit_cell(in.current_cell,input_data.evaluation_points[q]);

        const Quadrature<dim> quadrature_formula (quadrature_points);
        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values   |
                                 update_gradients |
                                 update_quadrature_points |
                                 update_JxW_values);

        fe_values.reinit(in.current_cell);
        this->get_material_model().fill_additional_material_model_inputs(in,
                                                                         this->get_solution(),
                                                                         fe_values,
                                                                         this->introspection());
        this->get_material_model().evaluate(in, out);

        if (this->get_parameters().formulation_temperature_equation
            == Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
          {
            // Overwrite the density by the reference density coming from the
            // adiabatic conditions as required by the formulation
            for (unsigned int q=0; q<n_quadrature_points; ++q)
              out.densities[q] = this->get_adiabatic_conditions().density(in.position[q]);
          }
        else if (this->get_parameters().formulation_temperature_equation
                 == Parameters<dim>::Formulation::TemperatureEquation::real_density)
          {
            // use real density
          }
        else
          AssertThrow(false, ExcNotImplemented());

        unsigned int index = 0;
        for (typename std::list<std::unique_ptr<HeatingModel::Interface<dim> > >::const_iterator
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
