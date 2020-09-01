/*
 Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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
          * @copydoc aspect::Particle::Property::Interface::initialize_one_particle_property()
          **/
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
          * @copydoc aspect::Particle::Property::Interface::update_particle_property()
          **/
          void
          update_particle_property (const unsigned int data_position,
                                    const Vector<double> &solution,
                                    const std::vector<Tensor<1,dim> > &gradients,
                                    typename ParticleHandler<dim>::particle_iterator &particle) const override;

          /**
          * @copydoc aspect::Particle::Property::Interface::need_update()
          **/
          UpdateTimeFlags
          need_update () const override;

          /**
          * @copydoc aspect::Particle::Property::Interface::get_needed_update_flags()
          **/
          UpdateFlags
          get_needed_update_flags () const override;

          /**
          * @copydoc aspect::Particle::Property::Interface::get_property_information()
          **/
          std::vector<std::pair<std::string, unsigned int> >
          get_property_information() const override;

        private:
          /**
           * Objects that are used to compute the particle property. Since the
           * object is expensive to create and is needed often it is kept as a
           * member variable. Because it is changed inside a const member function
           * (update_particle_property) it has to be mutable, but since it is
           * only used inside that function and always set before being used
           * that is not a problem. This implementation is not thread safe,
           * but it is currently not used in a threaded context.
           */
          mutable MaterialModel::MaterialModelInputs<dim> material_inputs;
          mutable MaterialModel::MaterialModelOutputs<dim> material_outputs;
      };
    }
  }
}

#endif

