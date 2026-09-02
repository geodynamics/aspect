/*
  Copyright (C) 2015 - 2026 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/particle/property/particle_generation_time.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      ParticleGenerationTime<dim>::initialize_one_particle_property(const Point<dim> &,
                                                                    std::vector<double> &data) const
      {
        const double current_time = this->get_time() / (this->convert_output_to_years() ? year_in_seconds : 1.0);
        data.push_back(current_time);
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      ParticleGenerationTime<dim>::get_property_information() const
      {
        return {{"particle generation time",1}};
      }



      template <int dim>
      InitializationModeForLateParticles
      ParticleGenerationTime<dim>::late_initialization_mode () const
      {
        return initialize;
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(ParticleGenerationTime,
                                        "particle generation time",
                                        "A plugin that stores the model time at which the particle was generated "
                                        "in the model domain.")
    }
  }
}
