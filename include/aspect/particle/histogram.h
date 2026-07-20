/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#ifndef _aspect_particle_histogram_h
#define _aspect_particle_histogram_h

#include <aspect/simulator_access.h>
#include <deal.II/particles/property_pool.h>
#include <deal.II/particles/particle_handler.h>


namespace aspect
{
  namespace Particle
  {
    /**
     * This class computes a spatial histogram describing the distribution of
     * particles in a cell. The number of spatial bins is equal to granularity^dimensions.
     */
    template <int dim>
    class ParticleHistogram
    {
      public:
        const unsigned int granularity;

        /**
         * This is the constructor for ParticleHistogram. ParticleHistogram sorts
         * the particles in a cell into different histogram bins depending on the
         * position of the particle in the cell. These histogram bins are compared to
         * a virtual histogram representing a uniform distribution of particles in the
         * cell and the difference between the real and virtual histograms is reported
         * as a score representing the degree of particle clustering in the cell.
         *
         * @param granularity determines the number spatial bins into which ParticleHistogram
         * divides each cell. The real number of bins is equal to granularity^dim, where dim
         * is the number of spatial dimensions, either 2 or 3.
         */
        ParticleHistogram(const unsigned int granularity);

        /**
         * Calculates and returns a score between 0 and 1 representing the level of clustering of the
         * particle range within the cell, where 0 represents a completely even distribution
         * and 1 represents the most clustered distribution that the histogram is able to
         * measure. The precision of the score depends on the ParticleHistogram's granularity
         * parameter.
         * @param particle_range is the particle range for which to evaluate the clustering
         * score
         * @param n_particles_in_cell is the number of particles in the particle range
         * @return A double ranging from 0 to 1 representing the degree to which the particles in
         * the provided particle range are clustered together, with 1 being the most clustered
         */
        double evaluate_particle_range(const typename Particles::ParticleHandler<dim>::particle_iterator_range &particle_range,
                                       const unsigned int n_particles_in_cell);

        /**
         * Sorts all of the particles within the cell into a deal.II table based on their position.
         * @param cell The cell for which to compute the particle distribution.
         * @param bucket_width The size (relative to the size of the cell) of each bucket in the table.
         * @return The table with the particle information.
         */
        Table<dim,unsigned int>
        sort_particles_into_buckets(const typename Particles::ParticleHandler<dim>::particle_iterator_range &particle_range,
                                    const double bucket_width) const;
    };
  }
}

#endif
