/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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

#include <aspect/particle/particle_handler.h>

namespace aspect
{
  namespace Particle
  {
    template <int dim,int spacedim>
    ParticleHandler<dim,spacedim>::ParticleHandler()
      :
      triangulation(),
      mpi_communicator(),
      global_number_of_particles(0),
      global_max_particles_per_cell(0),
      next_free_particle_index(0)
    {}



    template <int dim,int spacedim>
    ParticleHandler<dim,spacedim>::ParticleHandler(const parallel::distributed::Triangulation<dim,spacedim> &triangulation,
                                                   const MPI_Comm mpi_communicator)
      :
      triangulation(&triangulation, typeid(*this).name()),
      mpi_communicator(mpi_communicator),
      global_number_of_particles(0),
      global_max_particles_per_cell(0),
      next_free_particle_index(0)
    {}



    template <int dim,int spacedim>
    ParticleHandler<dim,spacedim>::~ParticleHandler()
    {}



    template <int dim,int spacedim>
    void
    ParticleHandler<dim,spacedim>::initialize(const parallel::distributed::Triangulation<dim,spacedim> &tria,
                                              const MPI_Comm communicator)
    {
      triangulation = &tria;
      mpi_communicator = communicator;
    }



    template <int dim,int spacedim>
    void
    ParticleHandler<dim,spacedim>::clear()
    {
      particles.clear();
      global_number_of_particles = 0;
      next_free_particle_index = 0;
      global_max_particles_per_cell = 0;
    }



    template <int dim,int spacedim>
    typename ParticleHandler<dim,spacedim>::particle_iterator
    ParticleHandler<dim,spacedim>::begin() const
    {
      return particle_iterator(particles,(const_cast<ParticleHandler<dim,spacedim> *> (this))->particles.begin());
    }



    template <int dim,int spacedim>
    typename ParticleHandler<dim,spacedim>::particle_iterator
    ParticleHandler<dim,spacedim>::begin()
    {
      return ParticleHandler<dim,spacedim>::particle_iterator(particles,particles.begin());
    }



    template <int dim,int spacedim>
    typename ParticleHandler<dim,spacedim>::particle_iterator
    ParticleHandler<dim,spacedim>::end() const
    {
      return (const_cast<ParticleHandler<dim,spacedim> *> (this))->end();
    }



    template <int dim,int spacedim>
    typename ParticleHandler<dim,spacedim>::particle_iterator
    ParticleHandler<dim,spacedim>::end()
    {
      return ParticleHandler<dim,spacedim>::particle_iterator(particles,particles.end());
    }



    template <int dim,int spacedim>
    boost::iterator_range<typename ParticleHandler<dim,spacedim>::particle_iterator>
    ParticleHandler<dim,spacedim>::particle_range_in_cell(const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell)
    {
      const types::LevelInd level_index = std::make_pair<int, int> (cell->level(),cell->index());
      const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::iterator,
            typename std::multimap<types::LevelInd, Particle<dim> >::iterator> particles_in_cell
            = particles.equal_range(level_index);

      return boost::make_iterator_range(particle_iterator(particles,particles_in_cell.first),
                                        particle_iterator(particles,particles_in_cell.second));
    }



    template <int dim,int spacedim>
    void
    ParticleHandler<dim,spacedim>::remove_particle(const ParticleHandler<dim,spacedim>::particle_iterator &particle)
    {
      particles.erase(particle->particle);
    }



    template <int dim,int spacedim>
    std::multimap<types::LevelInd, Particle<dim,spacedim> > &
    ParticleHandler<dim,spacedim>::get_particles()
    {
      return particles;
    }



    template <int dim,int spacedim>
    const std::multimap<types::LevelInd, Particle<dim,spacedim> > &
    ParticleHandler<dim,spacedim>::get_particles() const
    {
      return particles;
    }



    template <int dim,int spacedim>
    types::particle_index
    ParticleHandler<dim,spacedim>::n_global_particles() const
    {
      return global_number_of_particles;
    }



    template <int dim,int spacedim>
    types::particle_index
    ParticleHandler<dim,spacedim>::n_locally_owned_particles() const
    {
      return particles.size();
    }



    template <int dim,int spacedim>
    void
    ParticleHandler<dim,spacedim>::update_n_global_particles()
    {
      global_number_of_particles = dealii::Utilities::MPI::sum (particles.size(), mpi_communicator);
    }



    template <int dim,int spacedim>
    unsigned int
    ParticleHandler<dim,spacedim>::n_particles_in_cell(const typename Triangulation<dim>::active_cell_iterator &cell) const
    {
      const types::LevelInd found_cell = std::make_pair<int, int> (cell->level(),cell->index());
      return particles.count(found_cell);
    }



    template <int dim,int spacedim>
    void
    ParticleHandler<dim,spacedim>::update_next_free_particle_index()
    {
      types::particle_index locally_highest_index = 0;
      for (particle_iterator particle = begin(); particle != end(); ++particle)
        locally_highest_index = std::max(locally_highest_index,particle->get_id());

      next_free_particle_index = dealii::Utilities::MPI::max (locally_highest_index, mpi_communicator) + 1;
    }



    template <int dim,int spacedim>
    void
    ParticleHandler<dim,spacedim>::update_global_max_particles_per_cell()
    {
      unsigned int local_max_particles_per_cell(0);
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation->begin_active();
      for (; cell!=triangulation->end(); ++cell)
        if (cell->is_locally_owned())
          {
            const unsigned int particles_in_cell = n_particles_in_cell(cell);
            local_max_particles_per_cell = std::max(local_max_particles_per_cell,
                                                    particles_in_cell);
          }

      global_max_particles_per_cell = dealii::Utilities::MPI::max(local_max_particles_per_cell,mpi_communicator);
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class ParticleHandler<dim,dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
