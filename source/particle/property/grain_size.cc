/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/property/grain_size.h>
#include <aspect/material_model/grain_size.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      GrainSize<dim>::GrainSize ()
        :
        material_inputs(1,0),
        material_outputs(1,0)
      {}



      template <int dim>
      void
      GrainSize<dim>::initialize ()
      {
        CitationInfo::add("grainsize");

        material_inputs  = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());
        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

        AssertThrow(this->introspection().compositional_name_exists("grain_size"),
                    ExcMessage("This particle property only makes sense if "
                               "there is a compositional field named 'grain_size'."));

        grain_size_index = this->introspection().compositional_index_for_name("grain_size");
      }



      template <int dim>
      void
      GrainSize<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                       std::vector<double> &data) const
      {
        // Set the initial composition to the initial grain size.
        data.push_back(this->get_initial_composition_manager().initial_composition(position,grain_size_index));
      }



      template <int dim>
      void
      GrainSize<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                 typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        material_inputs  = MaterialModel::MaterialModelInputs<dim>(inputs.solution.size(), this->n_compositional_fields());
        material_outputs = MaterialModel::MaterialModelOutputs<dim>(inputs.solution.size(), this->n_compositional_fields());
        material_inputs.requested_properties = MaterialModel::MaterialProperties::reaction_terms;
        material_inputs.current_cell = inputs.current_cell;

        unsigned int p = 0;
        for (auto particle: particles)
          {
            // Make sure all particles are in the same cell
            Assert(particle.get_surrounding_cell() == inputs.current_cell,
                   ExcMessage("All particles must be in the same cell."));

            material_inputs.position[p] = particle.get_location();
            material_inputs.temperature[p] = inputs.solution[p][this->introspection().component_indices.temperature];
            material_inputs.pressure[p] = inputs.solution[p][this->introspection().component_indices.pressure];

            for (unsigned int d = 0; d < dim; ++d)
              material_inputs.velocity[p][d] = inputs.solution[p][this->introspection().component_indices.velocities[d]];

            for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
              material_inputs.composition[p][n] = inputs.solution[p][this->introspection().component_indices.compositional_fields[n]];

            material_inputs.composition[p][grain_size_index] = particle.get_properties()[this->data_position];

            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = inputs.gradients[p][this->introspection().component_indices.velocities[d]];

            material_inputs.strain_rate[p] = symmetrize (grad_u);

            ++p;
          }

        this->get_material_model().evaluate(material_inputs,
                                            material_outputs);

        p = 0;
        for (auto &particle: particles)
          {
            particle.get_properties()[this->data_position] += material_outputs.reaction_terms[p][grain_size_index];
            ++p;
          }
      }



      template <int dim>
      InitializationModeForLateParticles
      GrainSize<dim>::late_initialization_mode () const
      {
        return interpolate_respect_boundary;
      }



      template <int dim>
      UpdateTimeFlags
      GrainSize<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      GrainSize<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.velocities[component] == true)
          return update_values | update_gradients;

        return update_values;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      GrainSize<dim>::get_property_information() const
      {
        return {{"grain_size",1}};
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(GrainSize,
                                        "grain size",
                                        "A plugin in which the particle property is "
                                        "defined as the evolving grain size of a particle. "
                                        "See the grain_size material model "
                                        "documentation for more detailed information.")

    }
  }
}
