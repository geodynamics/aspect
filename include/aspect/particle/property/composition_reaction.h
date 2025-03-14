/*
 Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_property_composition_reaction_h
#define _aspect_particle_property_composition_reaction_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that initializes particle properties based on
       * the initial value of the compositional fields in
       * the model and updates them based on a description of
       * reactions between the fields.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class CompositionReaction : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Initialization function. This function is called once at the
           * creation of every particle for every property to initialize its
           * value.
           *
           * @param [in] position The current particle position.
           * @param [in,out] particle_properties The properties of the particle
           * that is initialized within the call of this function. The purpose
           * of this function should be to extend this vector by a number of
           * properties.
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
           * Returns an enum, which determines how this particle property is
           * initialized for particles that are created later than the initial
           * particle generation. For this property the value of
           * generated particles is interpolated from existing particles, unless
           * the particle is in a boundary cell that has a Dirichlet boundary
           * condition, in which case it uses the boundary condition value.
           */
          InitializationModeForLateParticles
          late_initialization_mode () const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::need_update()
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
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
           * A function object representing the area where each
           * reaction occurs.
           */
          std::unique_ptr<Functions::ParsedFunction<dim>> reaction_area;

          /**
           * The coordinate representation to evaluate the reaction_area
           * function. Possible choices are depth, cartesian and spherical.
           */
          Utilities::Coordinates::CoordinateSystem coordinate_system;

          /**
           * A function object representing the change in composition
           * during a reaction.
           * Inputs to this function are the reactant and the product
           * of the reaction, so that the values of the compositions
           * themselves can affect the reaction.
           */
          std::unique_ptr<Functions::ParsedFunction<2>> reaction_rate;

          /**
           * Vectors that store the indices of the compositional fields
           * taking part in each reaction.
           */
          std::vector<unsigned int> reactant_indices;
          std::vector<unsigned int> product_indices;

          /**
           * Vector of the times when each reaction should occur.
           */
          std::vector<double> reaction_times;

          /**
           * Reactions can include a 'background' field, meaning that
           * reactants vanish or products are created out of nowhere.
           * To facilitate this functionality, this (non-existing)
           * background field needs an index.
           */
          constexpr static unsigned int background_index = numbers::invalid_unsigned_int;
      };
    }
  }
}

#endif
