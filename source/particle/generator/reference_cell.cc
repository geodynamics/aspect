// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/particle/generator/reference_cell.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      ReferenceCell<dim>::generate_particles(Particles::ParticleHandler<dim> &particle_handler)
      {
        const std::vector<Point<dim>> reference_locations = generate_particle_positions_in_unit_cell();

        Particles::Generators::regular_reference_locations(this->get_triangulation(),
                                                           reference_locations,
                                                           particle_handler,
                                                           this->get_mapping());
      }


      template <int dim>
      std::vector<Point<dim>>
      ReferenceCell<dim>::generate_particle_positions_in_unit_cell()
      {
        std::vector<Point<dim>> particle_positions;
        std::array<double, dim> spacing;

        // Calculate separation of particles
        for (unsigned int i = 0; i < dim; ++i)
          spacing[i] = 1.0 / number_of_particles[i];

        for (unsigned int i = 0; i < number_of_particles[0]; ++i)
          {
            for (unsigned int j = 0; j < number_of_particles[1]; ++j)
              {
                if (dim == 2)
                  {
                    const Point<dim> position_unit = Point<dim>(i * spacing[0] + spacing[0] / 2,
                                                                j * spacing[1] + spacing[1] / 2);
                    particle_positions.push_back(position_unit);
                  }
                else if (dim == 3)
                  {
                    for (unsigned int k = 0; k < number_of_particles[2]; ++k)
                      {
                        const Point<dim> position_unit = Point<dim>(i * spacing[0] + spacing[0] / 2,
                                                                    j * spacing[1] + spacing[1] / 2,
                                                                    k * spacing[2] + spacing[2] / 2);
                        particle_positions.push_back(position_unit);
                      }
                  }
                else
                  ExcNotImplemented();
              }
          }

        return particle_positions;
      }


      template <int dim>
      void
      ReferenceCell<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Reference cell");
          {
            prm.declare_entry ("Number of particles per cell per direction", "2",
                               Patterns::List(Patterns::Integer(1)),
                               "List of number of particles to create per cell and spatial dimension. "
                               "The size of the list is the number of spatial dimensions. If only "
                               "one value is given, then each spatial dimension is set to the same value. "
                               "The list of numbers are parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      ReferenceCell<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Reference cell");
          {
            const auto n_particles_per_direction = Utilities::possibly_extend_from_1_to_N (
                                                     Utilities::string_to_int(
                                                       Utilities::split_string_list(prm.get("Number of particles per cell per direction"))),
                                                     dim,
                                                     "Number of particles per cell per direction");

            for (const auto &n_particle_direction: n_particles_per_direction)
              number_of_particles.push_back(static_cast<unsigned int> (n_particle_direction));
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(ReferenceCell,
                                         "reference cell",
                                         "Generates a uniform distribution of particles per cell and spatial direction in "
                                         "the unit cell and transforms each of the particles back to real region in the model "
                                         "domain. Uniform here means the particles will be generated with an equal spacing in "
                                         "each spatial dimension.")
    }
  }
}
