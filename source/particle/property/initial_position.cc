// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/global.h>
#include <aspect/particle/property/initial_position.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      InitialPosition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                             std::vector<double> &data) const
      {
        for (unsigned int i = 0; i < dim; ++i)
          data.push_back(position[i]);
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      InitialPosition<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("initial position",static_cast<unsigned int> (dim)));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(InitialPosition,
                                        "initial position",
                                        "Implementation of a plugin in which the particle "
                                        "property is given as the initial position "
                                        "of the particle. This property is vector-valued with "
                                        "as many components as there are space dimensions. "
                                        "In practice, it is often most useful to only "
                                        "visualize one of the components of this vector, "
                                        "or the magnitude of the vector. For example, in "
                                        "a spherical mantle simulation, the magnitude of this "
                                        "property equals the starting radius of a particle, "
                                        "and is thereby indicative of which part of the "
                                        "mantle a particle comes from.")
    }
  }
}
