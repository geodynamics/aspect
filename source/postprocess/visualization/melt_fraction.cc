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


#include <aspect/postprocess/visualization/melt_fraction.h>
#include <aspect/melt.h>
#include <aspect/material_model/reaction_model/pyroxenite_melting.h>
#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>

#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      MeltFraction<dim>::
      MeltFraction ()
        :
        DataPostprocessorScalar<dim> ("melt_fraction",
                                      update_values | update_quadrature_points),
        Interface<dim>("")  // a fraction, so physical units "1"
      {}



      template <int dim>
      void
      MeltFraction<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        // in case the material model computes the melt fraction itself, we use that output
        if (MaterialModel::MeltFractionModel<dim>::is_melt_fraction_model(this->get_material_model()))
          {
            MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                       this->introspection());
            MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                         this->n_compositional_fields());
            MeltHandler<dim>::create_material_model_outputs(out);

            // Compute the melt fraction...
            this->get_material_model().evaluate(in, out);

            std::vector<double> melt_fractions(n_quadrature_points);
            MaterialModel::MeltFractionModel<dim>::as_melt_fraction_model(this->get_material_model())
            .melt_fractions(in, melt_fractions, &out);

            for (unsigned int q=0; q<n_quadrature_points; ++q)
              computed_quantities[q](0) = melt_fractions[q];
          }
        else
          for (unsigned int q=0; q<n_quadrature_points; ++q)
            {
              const double pressure    = input_data.solution_values[q][this->introspection().component_indices.pressure];
              const double temperature = input_data.solution_values[q][this->introspection().component_indices.temperature];
              std::vector<double> composition(this->n_compositional_fields());

              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                composition[c] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[c]];

              // melting model 1: anhydrous melting of peridotite after Katz, 2003
              double peridotite_melt_fraction = katz2003_model.melt_fraction(temperature, pressure);

              // melting model 2: melting of pyroxenite after Sobolev et al., 2011 here. 
              // Can be replaced by any other melting model that is implemented in ASPECT with a "melt_fraction" function.
              double comp_2_melt_fraction = pyroxenite_model.melt_fraction(temperature, pressure);
              double melt_fraction;
              if (this->introspection().compositional_name_exists("composition_2"))
                {
                  const unsigned int comp_2_index = this->introspection().compositional_index_for_name("composition_2");
                  melt_fraction = composition[comp_2_index] * comp_2_melt_fraction +
                                  (1-composition[comp_2_index]) * peridotite_melt_fraction;
                }
              else
                melt_fraction = peridotite_melt_fraction;

              computed_quantities[q](0) = melt_fraction;
            }
      }



      template <int dim>
      void
      MeltFraction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt fraction");
            {
              MaterialModel::ReactionModel::Katz2003MantleMelting<dim>::declare_parameters(prm);
              MaterialModel::ReactionModel::PyroxeniteMelting<dim>::declare_parameters(prm);
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      MeltFraction<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt fraction");
            {
              // initialize simulator access of the katz2003_mantle_melting and pyroxenite_melting model
              katz2003_model.initialize_simulator(this->get_simulator());
              pyroxenite_model.initialize_simulator(this->get_simulator());
              
              // call parse_parameters of the katz2003_mantle_melting model and pyroxenite_melting to get the solidus and melt fraction parameters
              katz2003_model.parse_parameters(prm);
              pyroxenite_model.parse_parameters(prm);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltFraction,
                                                  "melt fraction", // TODO write down equations here
                                                  "A visualization output object that generates output "
                                                  "for the melt fraction at the temperature and "
                                                  "pressure of the current point. If the material model computes "
                                                  "a melt fraction, this is the quantity that will be visualized. "
                                                  "Otherwise, a specific parametrization for batch melting "
                                                  "(as described in the following) will be used. "
                                                  "It does not take into account latent heat. "
                                                  "If there are no compositional fields, or no fields called 'composition_2', "
                                                  " this postprocessor will visualize the melt fraction of peridotite "
                                                  "(calculated using the anhydrous model of Katz, 2003). "
                                                  "If there is a compositional field called 'composition_2', the "
                                                  "postprocessor assumes that this compositional "
                                                  "field is the content of pyroxenite, and will visualize "
                                                  "the melt fraction for a mixture of peridotite and pyroxenite "
                                                  "(using the melting model of Sobolev, 2011 for pyroxenite). "
                                                  "All the parameters that were used in these calculations "
                                                  "can be changed in the input file, the most relevant maybe "
                                                  "being the mass fraction of Cpx in peridotite in the Katz "
                                                  "melting model (Mass fraction cpx), which right now has a "
                                                  "default of 15\\%."
                                                  "The corresponding $p$-$T$-diagrams can be generated by running "
                                                  "the tests melt\\_postprocessor\\_peridotite and "
                                                  "melt\\_postprocessor\\_pyroxenite."
                                                  "\n\n"
                                                  "Physical units: None.")
    }
  }
}
