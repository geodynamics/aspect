/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/melt.h>
#include <aspect/melt.h>
#include <aspect/utilities.h>
#include <aspect/simulator.h>
#include <aspect/material_model/interface.h>

#include <deal.II/numerics/data_out.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
      MeltMaterialProperties<dim>::
      MeltMaterialProperties ()
        :
        DataPostprocessor<dim> (),
        // What is being output depends on run-time parameters. We can not
        // describe physical units to the base class at this point.
        Interface<dim>("")
      {}



      template <int dim>
      std::vector<std::string>
      MeltMaterialProperties<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;

        for (const auto &property_name : property_names)
          if (property_name == "fluid density gradient")
            for (unsigned int i=0; i<dim; ++i)
              solution_names.emplace_back("fluid_density_gradient");
          else if (property_name == "compaction pressure")
            {
              solution_names.emplace_back("p_c");
            }
          else
            {
              solution_names.push_back(property_name);
              std::replace(solution_names.back().begin(),solution_names.back().end(),' ', '_');
            }

        return solution_names;
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      MeltMaterialProperties<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
        for (const auto &property_name : property_names)
          {
            if (property_name == "fluid density gradient")
              {
                for (unsigned int c=0; c<dim; ++c)
                  interpretation.push_back (DataComponentInterpretation::component_is_part_of_vector);
              }
            else
              interpretation.push_back (DataComponentInterpretation::component_is_scalar);
          }

        return interpretation;
      }



      template <int dim>
      UpdateFlags
      MeltMaterialProperties<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_quadrature_points;
      }



      template <int dim>
      void
      MeltMaterialProperties<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true,
                    ExcMessage("'Include melt transport' has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());

        // Set use_strain_rates to true since the compaction viscosity might also depend on the strain rate.
        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        MeltHandler<dim>::create_material_model_outputs(out);

        this->get_material_model().evaluate(in, out);
        MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim>>();
        AssertThrow(melt_outputs != nullptr,
                    ExcMessage("Need MeltOutputs from the material model for computing the melt properties."));

        const double p_c_scale = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>(this->get_material_model()).p_c_scale(in,
                                 out,
                                 this->get_melt_handler(),
                                 true);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            unsigned output_index = 0;
            for (unsigned int i=0; i<property_names.size(); ++i, ++output_index)
              {
                if (property_names[i] == "compaction viscosity")
                  computed_quantities[q][output_index] = melt_outputs->compaction_viscosities[q];
                else if (property_names[i] == "fluid viscosity")
                  computed_quantities[q][output_index] = melt_outputs->fluid_viscosities[q];
                else if (property_names[i] == "permeability")
                  computed_quantities[q][output_index] = melt_outputs->permeabilities[q];
                else if (property_names[i] == "fluid density")
                  computed_quantities[q][output_index] = melt_outputs->fluid_densities[q];
                else if (property_names[i] == "fluid density gradient")
                  {
                    for (unsigned int k=0; k<dim; ++k, ++output_index)
                      {
                        computed_quantities[q][output_index] = melt_outputs->fluid_density_gradients[q][k];
                      }
                    --output_index;
                  }
                else if (property_names[i] == "compaction pressure")
                  {
                    const unsigned int pc_comp_idx = this->introspection().variable("compaction pressure").first_component_index;
                    const double p_c_bar = input_data.solution_values[q][pc_comp_idx];

                    computed_quantities[q][output_index] = p_c_scale * p_c_bar;
                  }
                else if (property_names[i] == "darcy coefficient")
                  {
                    const double K_D = this->get_melt_handler().limited_darcy_coefficient(melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q], p_c_scale > 0);
                    computed_quantities[q][output_index] = K_D;
                  }
                else if (property_names[i] == "darcy coefficient no cutoff")
                  {
                    const double K_D_no_cut = melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q];
                    computed_quantities[q][output_index] = K_D_no_cut;
                  }
                else if (property_names[i] == "is melt cell")
                  {
                    computed_quantities[q][output_index] = this->get_melt_handler().is_melt_cell(in.current_cell)? 1.0 : 0.0;
                  }
                else if (property_names[i] == "compaction length")
                  {
                    const double compaction_length = std::sqrt((out.viscosities[q] + 4./3. * melt_outputs->compaction_viscosities[q])
                                                               * melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]);
                    computed_quantities[q][output_index] = compaction_length;
                  }
                else
                  AssertThrow(false, ExcNotImplemented());
              }
          }
      }



      template <int dim>
      void
      MeltMaterialProperties<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt material properties");
            {
              const std::string pattern_of_names
                = "compaction viscosity|fluid viscosity|permeability|"
                  "fluid density|fluid density gradient|is melt cell|"
                  "darcy coefficient|darcy coefficient no cutoff|"
                  "compaction length";

              prm.declare_entry("List of properties",
                                "compaction viscosity,permeability",
                                Patterns::MultipleSelection(pattern_of_names),
                                "A comma separated list of melt properties that should be "
                                "written whenever writing graphical output. "
                                "The following material properties are available:\n\n"
                                +
                                pattern_of_names);
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      MeltMaterialProperties<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt material properties");
            {
              property_names = Utilities::split_string_list(prm.get ("List of properties"));
              AssertThrow(Utilities::has_unique_entries(property_names),
                          ExcMessage("The list of strings for the parameter "
                                     "'Postprocess/Visualization/Melt material properties/List of properties' contains entries more than once. "
                                     "This is not allowed. Please check your parameter file."));

              // Always output compaction pressure
              property_names.insert(property_names.begin(),"compaction pressure");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltMaterialProperties,
                                                  "melt material properties",
                                                  "A visualization output object that generates output "
                                                  "for melt related properties of the material model. Note "
                                                  "that this postprocessor always outputs the compaction pressure, "
                                                  "but can output a large range of additional properties, as "
                                                  "selected in the ``List of properties'' parameter."
                                                  "\n\n"
                                                  "Physical units: Various, depending on what is being output.")

    }
  }
}
