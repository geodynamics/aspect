/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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
      ReferenceCell<dim>::generate_particles(std::multimap<Particles::internal::LevelInd, Particle<dim> > &particles)
      {
        const std::vector<Point<dim> > particles_in_unit_cell = generate_particle_positions_in_unit_cell();

        types::particle_index n_particles_to_generate = this->get_triangulation().n_locally_owned_active_cells() * particles_in_unit_cell.size();
        types::particle_index prefix_sum = 0;

#if DEAL_II_VERSION_GTE(9,1,0)
        const int ierr = MPI_Scan(&n_particles_to_generate, &prefix_sum, 1, DEAL_II_PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());
#else
        const int ierr = MPI_Scan(&n_particles_to_generate, &prefix_sum, 1, PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());
#endif
        AssertThrowMPI(ierr);

        types::particle_index particle_index = prefix_sum - n_particles_to_generate;

        for (const auto &cell : this->get_triangulation().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              for (const auto &itr_particles_in_unit_cell : particles_in_unit_cell)
                {
                  const Point<dim> position_real = this->get_mapping().transform_unit_to_real_cell(cell,
                                                   itr_particles_in_unit_cell);
                  const Particle<dim> particle(position_real, itr_particles_in_unit_cell, particle_index);
                  const Particles::internal::LevelInd cell_index(cell->level(), cell->index());
                  particles.insert(std::make_pair(cell_index, particle));
                  ++particle_index;
                }
            }
      }


      template <int dim>
      std::vector<Point<dim> >
      ReferenceCell<dim>::generate_particle_positions_in_unit_cell()
      {
        std::vector<Point<dim> > particle_positions;
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
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell");
              {
                prm.declare_entry ("Number of particles per cell per direction", "2",
                                   Patterns::List(Patterns::Double (0.)),
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
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      ReferenceCell<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell");
              {
                std::vector<double> n_particles_tmp = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Number of particles per cell per direction"))),
                                                                                              dim,
                                                                                              "Number of particles per cell per direction");
                for (double itr : n_particles_tmp)
                  number_of_particles.push_back(static_cast<unsigned int>(itr));
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
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
                                         "each spatial dimension")
    }
  }
}
