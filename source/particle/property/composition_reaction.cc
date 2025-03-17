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

#include <aspect/particle/property/composition_reaction.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      CompositionReaction<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                                 std::vector<double> &data) const
      {
        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          data.push_back(this->get_initial_composition_manager().initial_composition(position,i));
      }



      template <int dim>
      void
      CompositionReaction<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &/*inputs*/,
                                                           typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        // we get time passed as seconds (always) but may want
        // to reinterpret it in years
        if (this->convert_output_to_years())
          reaction_rate->set_time (this->get_time() / year_in_seconds);
        else
          reaction_rate->set_time (this->get_time());

        // Check if any reaction occurs during the current time step.
        const unsigned int n_reactions = reactant_indices.size();
        std::vector<bool> reaction_occurs (n_reactions, false);

        for (unsigned int i=0; i<n_reactions; ++i)
          if (this->get_time() >= reaction_times[i] &&
              (this->get_time() - this->get_timestep()) < reaction_times[i])
            reaction_occurs[i] = true;

        if (std::find(reaction_occurs.begin(), reaction_occurs.end(), true) == reaction_occurs.end())
          return;

        // Loop over all particles to apply reactions.
        for (auto &particle: particles)
          {
            const Utilities::NaturalCoordinate<dim> point =
              this->get_geometry_model().cartesian_to_other_coordinates(particle.get_location(), coordinate_system);
            const ArrayView<double> particle_properties = particle.get_properties();

            // Loop over all reactions and compute total change for each compositional field.
            for (unsigned int i=0; i<n_reactions; ++i)
              {
                std::array<double,2> reactant_and_product;
                if (reactant_indices[i] == background_index)
                  reactant_and_product[0] = 0.0;
                else
                  reactant_and_product[0] = particle.get_properties()[this->data_position + reactant_indices[i]];

                if (product_indices[i] == background_index)
                  reactant_and_product[1] = 0.0;
                else
                  reactant_and_product[1] = particle.get_properties()[this->data_position + product_indices[i]];

                const double area = reaction_area->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), i);
                const double rate = reaction_rate->value(Utilities::convert_array_to_point<2>(reactant_and_product), i);

                // Subtract the reaction value from the reactant and add it to the product.
                // For the background, we do not have to do anything.
                if (reaction_occurs[i] && reactant_indices[i] != background_index)
                  particle_properties[this->data_position + reactant_indices[i]] -= area * rate;
                if (reaction_occurs[i] && product_indices[i] != background_index)
                  particle_properties[this->data_position + product_indices[i]] += area * rate;
              }
          }
      }



      template <int dim>
      InitializationModeForLateParticles
      CompositionReaction<dim>::late_initialization_mode () const
      {
        return interpolate_respect_boundary;
      }



      template <int dim>
      UpdateTimeFlags
      CompositionReaction<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      CompositionReaction<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;

        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          {
            std::ostringstream field_name;
            field_name << this->introspection().name_for_compositional_index(i) << " reaction";
            property_information.emplace_back(field_name.str(),1);
          }

        return property_information;
      }


      template <int dim>
      void
      CompositionReaction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Composition reaction");
        {
          prm.declare_entry ("List of reactants", "",
                             Patterns::List(Patterns::Anything()),
                             "Select the compositional fields that are the reaction inputs. "
                             "Each entry represents the input for one reaction. The parameter "
                             "should contain a list of compositional field names, one per reaction. "
                             "'background' can be selected to set up a reaction without reactants. "
                             "The length of this list determines the number of components in the "
                             "reaction function.");
          prm.declare_entry ("List of products", "",
                             Patterns::List(Patterns::Anything()),
                             "Select the compositional fields that are the reaction outputs. "
                             "Each entry represents the output of one reaction. The parameter "
                             "should contain a list of compositional field names, one per reaction. "
                             "'background' can be selected to set up a reaction without products. "
                             "Needs to have as many entries as the 'List of reactants'.");
          prm.declare_entry ("List of reaction times", "",
                             Patterns::List(Patterns::Double(0.)),
                             "List a specific point in time when each reaction should occur during "
                             "the simulation. If set to zero, the reaction occurs throughout the "
                             "whole simulation."
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");

          prm.enter_subsection("Reaction area function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);

            prm.declare_entry ("Coordinate system", "cartesian",
                               Patterns::Selection ("depth|cartesian|spherical"),
                               "A selection that determines the assumed coordinate "
                               "system for the function variables. Allowed values "
                               "are `depth', `cartesian' and `spherical'. `depth' "
                               "will create a function with only the first variable "
                               "being is non-zero, and this first variable is interpreted "
                               "as the depth of the point. `spherical' coordinates "
                               "are interpreted as r,phi or r,phi,theta in 2d/3d, "
                               "respectively, with theta being the polar angle.");
          }
          prm.leave_subsection();

          prm.enter_subsection("Reaction rate function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);

            prm.declare_entry ("Variable names", "reactant,product,t",
                               Patterns::Anything(),
                               "The names of the variables as they will be used in the function, "
                               "separated by commas. Instead of spatial coordinates, the inputs for "
                               "this function represent the value of the `reactant' and `product' "
                               "compositions so that the value of the current composition can be "
                               "used to calculate the change due to the reaction. Additionally, "
                               "`t' represents the time. You can use these variable "
                               "names in your function expression and they will be replaced by the "
                               "values of these variables at which the function is currently evaluated. "
                               "However, you can also choose a different set of names for the "
                               "independent variables at which to evaluate your function expression.");
            prm.declare_entry ("Function expression","1.",
                               Patterns::Anything(),
                               "Expression for the change in value of the reactant and reaction product "
                               "for each reaction, in dependence of the current values of reactant and "
                               "reaction product. One expression for each reaction should be listed, "
                               "separated by semicolons.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      CompositionReaction<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow(this->n_compositional_fields() > 0,
                    ExcMessage("You have requested the particle property <composition "
                               "reaction>, but the number of compositional fields is 0. "
                               "Please add compositional fields to your model, or remove "
                               "this particle property."));

        prm.enter_subsection("Composition reaction");
        {
          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
          compositional_field_names.insert(compositional_field_names.begin(),"background");

          std::vector<std::string> reactants = Utilities::split_string_list(prm.get("List of reactants"));
          std::vector<std::string> products = Utilities::split_string_list(prm.get("List of products"));
          reaction_times = Utilities::string_to_double
                           (Utilities::split_string_list(prm.get ("List of reaction times")));

          if (this->convert_output_to_years() == true)
            for (unsigned int i=0; i<reaction_times.size(); ++i)
              reaction_times[i] *= year_in_seconds;

          AssertThrow ((reactants.size() == products.size()),
                       ExcMessage ("The list of reactants and products need to have the "
                                   "same size. You provided "
                                   + Utilities::to_string(reactants.size()) + " reactants and "
                                   + Utilities::to_string(products.size()) + " products."));

          AssertThrow ((reactants.size() == reaction_times.size()),
                       ExcMessage ("There must be as many reaction times given "
                                   "as there are reactions. You provided "
                                   + Utilities::to_string(reaction_times.size()) +
                                   " reaction times but there are "
                                   + Utilities::to_string(reactants.size()) + " reactions."));

          reactant_indices.resize(reactants.size());
          product_indices.resize(products.size());

          for (unsigned int i=0; i<reactants.size(); ++i)
            {
              AssertThrow (std::find(compositional_field_names.begin(), compositional_field_names.end(), reactants[i])
                           != compositional_field_names.end(),
                           ExcMessage ("The reactant " + reactants[i]  + " you provided "
                                       "in the 'List of reactants' is not a compositional field name."));

              if (reactants[i] == "background")
                reactant_indices[i] = background_index;
              else
                reactant_indices[i] = this->introspection().compositional_index_for_name(reactants[i]);

              AssertThrow (std::find(compositional_field_names.begin(), compositional_field_names.end(), products[i])
                           != compositional_field_names.end(),
                           ExcMessage ("The product " + products[i]  + " you provided "
                                       "in the 'List of products' is not a compositional field name."));

              if (products[i] == "background")
                product_indices[i] = background_index;
              else
                product_indices[i] = this->introspection().compositional_index_for_name(products[i]);
            }

          prm.enter_subsection("Reaction area function");
          {

            coordinate_system = ::aspect::Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));

            try
              {
                reaction_area = std::make_unique<Functions::ParsedFunction<dim>>(reactants.size());
                reaction_area->parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Particles.CompositionReaction.ReactionAreaFunction'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'";
                throw;
              }
          }
          prm.leave_subsection();

          prm.enter_subsection("Reaction rate function");
          {
            try
              {
                reaction_rate = std::make_unique<Functions::ParsedFunction<2>>(reactants.size());
                reaction_rate->parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Particles.CompositionReaction.ReactionRateFunction'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'";
                throw;
              }
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
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(CompositionReaction,
                                        "composition reaction",
                                        "Implementation of a plugin in which the particle "
                                        "property is given as the initial composition "
                                        "at the particle's initial position, and is updated "
                                        "during the simulation time according to reactions "
                                        "that are specified as functions in the input file. "
                                        "Each reaction has exactly one reactant and one product. "
                                        "Each particle gets as many properties as there are "
                                        "compositional fields. "
                                        "The reactions are described by two functions, and "
                                        "the change in each composition at the time a reaction "
                                        "occurs is computed as the product of the two functions. "
                                        "The 'Reaction area function' describes the area "
                                        "where the reaction takes place. It can be spatially "
                                        "variable, but does not depend on time. The 'Reaction "
                                        "rate function' describes how the change in composition "
                                        "depends on these compositions themselves and on time. "
                                        "To use this particle property for a given compositional "
                                        "field, set the 'Mapped particle properties' to "
                                        "'name_of_field:name_of_field reaction', i.e., the "
                                        "name of the particle property for each field is "
                                        "the name of the compositional field with the word "
                                        "'reaction' added at the end.")
    }
  }
}
