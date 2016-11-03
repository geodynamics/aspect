/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/generator/reference_cell.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      ReferenceCell<dim>::ReferenceCell()
      {
        starting_particle_index = 0;
      }

      template <int dim>
      void
      ReferenceCell<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        const std::vector<Point<dim>> particles_in_unit_cell = generate_particle_positions_in_unit_cell();

        unsigned int n_particles_to_generate = this->get_triangulation().n_locally_owned_active_cells() * particles_in_unit_cell.size();
        unsigned int prefix_sum = 0;

        MPI_Scan(&n_particles_to_generate, &prefix_sum, 1, MPI_UNSIGNED, MPI_SUM, this->get_mpi_communicator());

        starting_particle_index += prefix_sum - n_particles_to_generate;
        types::particle_index particle_index = starting_particle_index;
        typename Triangulation<dim>::active_cell_iterator
        cell = this->get_triangulation().begin_active(),
        endc = this->get_triangulation().end();

        for (; cell != endc; cell++)
          {
            if (cell->is_locally_owned())
              {
                for (typename std::vector<Point<dim>>::const_iterator itr_particles_in_unit_cell = particles_in_unit_cell.begin();
                     itr_particles_in_unit_cell != particles_in_unit_cell.end();
                     itr_particles_in_unit_cell++)
                  {
                    try
                      {
                        Point<dim> position_real = this->get_mapping().transform_unit_to_real_cell(cell,
                                                                                                   *itr_particles_in_unit_cell);
                        const Particle<dim> particle(position_real, *itr_particles_in_unit_cell, particle_index);
                        const types::LevelInd cell_index(cell->level(), cell->index());
                        particles.insert(std::make_pair(cell_index, particle));
                        particle_index++;
                      }
                    catch (typename Mapping<dim>::ExcTransformationFailed &)
                      {
                        AssertThrow (true,
                                     ExcMessage("Couldn't generate particle (unusual cell shape?). "));
                      }
                  }
              }
          }

        starting_particle_index = Utilities::MPI::max(particle_index, this->get_mpi_communicator());
      }


      template <int dim>
      const std::vector<Point<dim>>
                                 ReferenceCell<dim>::generate_particle_positions_in_unit_cell()
      {
        std::vector<Point<dim>> particle_positions;
        std_cxx11::array<double, dim> spacing;

        // Calculate separation of particles
        for (unsigned int i = 0; i < dim; ++i)
          spacing[i] = 1.0 / number_of_particles[i];

        for (unsigned int i = 0; i < number_of_particles[0]; ++i)
          {
            for (unsigned int j = 0; j < number_of_particles[1]; ++j)
              {
                if (dim == 2)
                  {
                    Point<dim> position_unit = Point<dim>(i * spacing[0] + spacing[0] / 2,
                                                          j * spacing[1] + spacing[1] / 2);
                    particle_positions.push_back(position_unit);
                  }
                else if (dim == 3)
                  {
                    for (unsigned int k = 0; k < number_of_particles[2]; ++k)
                      {
                        Point<dim> position_unit = Point<dim>(i * spacing[0] + spacing[0] / 2,
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
          prm.enter_subsection("Tracers");
          {
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell");
              {
                prm.declare_entry ("Number of tracers per cell per direction", "10",
                                   Patterns::Double (0),
                                   "Total number of tracers to create for each cell for each spatial dimension."
                                   "The number is parsed as a floating point number (so that one can "
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
          prm.enter_subsection("Tracers");
          {
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell");
              {
                for (int i=0; i<dim; i++)
                  number_of_particles[i] = (unsigned int) (prm.get_double ("Number of tracers per cell per direction"));
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
                                         "Generates a uniform distribution of particles "
                                         "in the unit domain and transforms the locations back onto the real domain.")
    }
  }
}
