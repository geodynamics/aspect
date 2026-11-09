// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/particle/property/composition.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      Composition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                         std::vector<double> &data) const
      {
        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          data.push_back(this->get_initial_composition_manager().initial_composition(position,i));
      }



      template <int dim>
      void
      Composition<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                   typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        unsigned int p = 0;
        const auto &composition_components = this->introspection().component_indices.compositional_fields;
        for (auto &particle: particles)
          {
            for (unsigned int j = 0; j < this->n_compositional_fields(); ++j)
              {
                particle.get_properties()[this->data_position+j] = inputs.solution[p][composition_components[j]];
              }
            ++p;
          }
      }



      template <int dim>
      UpdateTimeFlags
      Composition<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      Composition<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.compositions[component] == true)
          return update_values;

        return update_default;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      Composition<dim>::get_property_information() const
      {

        AssertThrow(this->n_compositional_fields() > 0,
                    ExcMessage("You have requested the particle property <composition>, "
                               "but the number of compositional fields is 0. "
                               "Please add compositional fields to your model, or remove "
                               "this particle property."));

        std::vector<std::pair<std::string,unsigned int>> property_information;



        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          {
            const std::string field_name = this->introspection().name_for_compositional_index(i);
            property_information.emplace_back(field_name,1);
          }
        return property_information;
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(Composition,
                                        "composition",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined by the compositional fields in "
                                        "the model. This can be used to track solid composition "
                                        "evolution over time.")
    }
  }
}
