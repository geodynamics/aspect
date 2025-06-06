/*
 Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_property_elastic_stress_h
#define _aspect_particle_property_elastic_stress_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that calculates the viscoelastic stress a particle
       * has accumulated through time.
       * The implementation of this property is equivalent to the implementation
       * for compositional fields that is described in the viscoelastic material
       * model and elastic stress rheology module.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class ElasticStress : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ElasticStress ();

          void initialize () override;

          /**
           * Function to update particles after they have been restored
           * to their position and values from the beginning of the timestep.
           * This restoring happens at the beginning of nonlinear iterations
           * of iterative Advection solver schemes.
           */
          void
          update_particles (typename Particle::Manager<dim> &particle_manager) const;

          /**
           * @copydoc aspect::Particle::Property::Interface::initialize_one_particle_property()
           */
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::update_particle_properties()
           */
          void
          update_particle_properties (const ParticleUpdateInputs<dim> &inputs,
                                      typename ParticleHandler<dim>::particle_iterator_range &particles) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::need_update()
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::get_update_flags()
           */
          UpdateFlags
          get_update_flags (const unsigned int component) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::get_property_information()
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Objects that are used to compute the particle property. Since the
           * object is expensive to create and is needed often it is kept as a
           * member variable. Because it is changed inside a const member function
           * (update_particle_property) it has to be mutable, but since it is
           * only used inside that function and always set before being used
           * that is not a problem. This implementation is not thread safe,
           * but it is currently not used in a threaded context.
           *
           * The first two objects provide material model in- and outputs for
           * one particle at a time; the second set for all the particles
           * in a cell at the same time. TODO use one set for both?
           */
          mutable MaterialModel::MaterialModelInputs<dim> material_inputs;
          mutable MaterialModel::MaterialModelOutputs<dim> material_outputs;
          mutable MaterialModel::MaterialModelInputs<dim> material_inputs_cell;
          mutable MaterialModel::MaterialModelOutputs<dim> material_outputs_cell;

          /**
           * The indices of the compositional fields that represent components of the
           * viscoelastic stress tensors.
           */
          std::vector<unsigned int> stress_field_indices;

          /**
           * The indices of the compositional fields that do not represent components of the
           * viscoelastic stress tensors.
           */
          std::vector<unsigned int> non_stress_field_indices;

          /**
           * The weight given to the stress values stored on the particles in the
           * weighted average with the stress values interpolated from the compositional
           * fields to the particle location. The default value of 1 is more accurate,
           * but can be less stable.
           */
          double particle_weight;
      };
    }
  }
}

#endif
