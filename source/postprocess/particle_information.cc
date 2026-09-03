/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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


#include <aspect/particle/manager.h>
#include <aspect/postprocess/particle_information.h>

#include <limits>
#include <sstream>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string, std::string>
    ParticleInformation<dim>::execute(TableHandler &)
    {
      // Only report once:
      if (this->get_timestep_number() > 0 || this->get_pre_refinement_step() != std::numeric_limits<unsigned int>::max())
        return {"", ""};

      std::ostringstream output;
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
        {
          output << "Particle manager " << particle_manager_index + 1 << " properties:" << std::endl;

          const auto property_names =
            this->get_particle_manager(particle_manager_index).get_property_manager().get_data_info().get_property_names();
          for (const auto &property_name : property_names)
            // The particle integrator adds internal storage to every world.
            // This is not a particle property selected by the user.
            if (property_name != "internal: integrator properties")
              output << "  " << property_name << std::endl;
        }

      return {"Particle information:", output.str()};
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleInformation,
                                  "particle information",
                                  "A postprocessor that prints the particle properties in every "
                                  "particle manager at time zero.")
  }
}
