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
