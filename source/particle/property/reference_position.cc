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

#include <aspect/particle/property/reference_position.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      ReferencePosition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                               std::vector<double> &data) const
      {
        for (unsigned int i = 0; i < dim; ++i)
          data.push_back(position[i]);
      }

      template <int dim>
      void
      ReferencePosition<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &/*inputs*/,
                                                         typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        for (auto &particle: particles)
          for (unsigned int i = 0; i < dim; ++i)
            particle.get_properties()[this->data_position+i] = particle.get_reference_location()[i];
      }

      template <int dim>
      UpdateTimeFlags
      ReferencePosition<dim>::need_update() const
      {
        return update_output_step;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      ReferencePosition<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("reference position",dim));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(ReferencePosition,
                                        "reference position",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as the current reference position.")
    }
  }
}
