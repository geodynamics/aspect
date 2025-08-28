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

#include <aspect/particle/generator/uniform_box.h>

#include <array>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      UniformBox<dim>::generate_particles(Particles::ParticleHandler<dim> &particle_handler)
      {
        const Tensor<1,dim> P_diff = P_max - P_min;

        double volume(1.0);
        for (unsigned int i = 0; i < dim; ++i)
          volume *= P_diff[i];

        std::array<unsigned int,dim> n_particles_per_direction;
        std::array<double,dim> spacing;

        // Calculate separation of particles
        for (unsigned int i = 0; i < dim; ++i)
          {
#if DEAL_II_VERSION_GTE(9,6,0)
            n_particles_per_direction[i] = static_cast<unsigned int>(std::round(std::pow(n_particles * Utilities::pow(P_diff[i],dim) / volume, 1./dim)));
#else
            n_particles_per_direction[i] = static_cast<unsigned int>(std::round(std::pow(n_particles * std::pow(P_diff[i],dim) / volume, 1./dim)));
#endif

            spacing[i] = P_diff[i] / fmax(n_particles_per_direction[i] - 1,1);
          }

        types::particle_index particle_index = 0;

        for (unsigned int i = 0; i < n_particles_per_direction[0]; ++i)
          {
            for (unsigned int j = 0; j < n_particles_per_direction[1]; ++j)
              {
                if (dim == 2)
                  {
                    const Point<dim> particle_position = Point<dim> (P_min[0]+i*spacing[0],P_min[1]+j*spacing[1]);
                    this->insert_particle_at_position(particle_position, particle_index, particle_handler);
                    ++particle_index;
                  }
                else if (dim == 3)
                  for (unsigned int k = 0; k < n_particles_per_direction[2]; ++k)
                    {
                      const Point<dim> particle_position = Point<dim> (P_min[0]+i*spacing[0],P_min[1]+j*spacing[1],P_min[2]+k*spacing[2]);
                      this->insert_particle_at_position(particle_position, particle_index, particle_handler);
                      ++particle_index;
                    }

              }
          }

        particle_handler.update_cached_numbers();
      }


      template <int dim>
      void
      UniformBox<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Uniform box");
          {
            prm.declare_entry ("Number of particles", "1000",
                               Patterns::Double (0.),
                               "Total number of particles to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");

            prm.declare_entry ("Minimum x", "0.",
                               Patterns::Double (),
                               "Minimum x coordinate for the region of particles.");
            prm.declare_entry ("Maximum x", "1.",
                               Patterns::Double (),
                               "Maximum x coordinate for the region of particles.");
            prm.declare_entry ("Minimum y", "0.",
                               Patterns::Double (),
                               "Minimum y coordinate for the region of particles.");
            prm.declare_entry ("Maximum y", "1.",
                               Patterns::Double (),
                               "Maximum y coordinate for the region of particles.");
            prm.declare_entry ("Minimum z", "0.",
                               Patterns::Double (),
                               "Minimum z coordinate for the region of particles.");
            prm.declare_entry ("Maximum z", "1.",
                               Patterns::Double (),
                               "Maximum z coordinate for the region of particles.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      UniformBox<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Uniform box");
          {
            n_particles    = static_cast<types::particle_index>(prm.get_double ("Number of particles"));

            P_min(0) = prm.get_double ("Minimum x");
            P_max(0) = prm.get_double ("Maximum x");
            P_min(1) = prm.get_double ("Minimum y");
            P_max(1) = prm.get_double ("Maximum y");

            AssertThrow(P_min(0) < P_max(0), ExcMessage("Minimum x must be less than maximum x"));
            AssertThrow(P_min(1) < P_max(1), ExcMessage("Minimum y must be less than maximum y"));

            if (dim == 3)
              {
                P_min(2) = prm.get_double ("Minimum z");
                P_max(2) = prm.get_double ("Maximum z");

                AssertThrow(P_min(2) < P_max(2), ExcMessage("Minimum z must be less than maximum z"));
              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      ASPECT_REGISTER_PARTICLE_GENERATOR(UniformBox,
                                         "uniform box",
                                         "Generate a uniform distribution of particles "
                                         "over a rectangular domain in 2d or 3d. Uniform here means "
                                         "the particles will be generated with an equal spacing in "
                                         "each spatial dimension. Note that in order "
                                         "to produce a regular distribution the number of generated "
                                         "particles might not exactly match the one specified in the "
                                         "input file.")
    }
  }
}
