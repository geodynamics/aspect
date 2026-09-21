// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_property_viscoplastic_strain_invariant_h
#define _aspect_particle_property_viscoplastic_strain_invariant_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/visco_plastic.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that calculates the finite strain invariant a particle has
       * experienced and assigns it to either the plastic and/or viscous strain field based
       * on whether the material is plastically yielding, or the total strain field
       * used in the visco plastic material model. The implementation of this property
       * is equivalent to the implementation for compositional fields that is located in
       * the plugin <code>benchmarks/buiter_et_al_2008_jgr/plugin/finite_strain_invariant.cc</code>,
       * and is effectively the same as what the visco plastic material model uses for compositional fields.
       * @ingroup ParticleProperties
       */
      template <int dim>
      class ViscoPlasticStrainInvariant : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ViscoPlasticStrainInvariant ();

          void initialize () override;


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

        private:
          unsigned int n_components;

          /**
           * An object that is used to compute the particle property. Since the
           * object is expensive to create and is needed often it is kept as a
           * member variable. Because it is changed inside a const member function
           * (update_particle_property) it has to be mutable, but since it is
           * only used inside that function and always set before being used
           * that is not a problem. This implementation is not thread safe,
           * but it is currently not used in a threaded context.
           */
          mutable MaterialModel::MaterialModelInputs<dim> material_inputs;
      };
    }
  }
}

#endif
