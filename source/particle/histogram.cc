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

#include <aspect/particle/histogram.h>
#include <aspect/particle/manager.h>
#include <deal.II/base/table.h>

namespace aspect
{
  namespace Particle
  {

    template <int dim>
    ParticleHistogram<dim>::ParticleHistogram(const unsigned int granularity)
      :
      granularity (granularity)
    {}



    template <int dim>
    double ParticleHistogram<dim>::evaluate_particle_range(const typename Particles::ParticleHandler<dim>::particle_iterator_range &particle_range,
                                                           const unsigned int n_particles_in_cell)
    {
      /*
        These are the steps taken to compute the density score:
        Create a table with perfect distribution containing an even number of particles in each bucket.
        Create a table with the "worst case" distribution in which all particles are in one bucket.
        Take the distance squared between the "perfect" table and the worst case table, treating both tables as vectors with granularity^dim elements.
        Create a table with the actual distribution.
        Take distance between "perfect" table and the actual table.
        Return the ratio between the observed distance and the worst case distance.
        A value of 1 represents the worst possible distribution, a value of 0 represents a completely even distribution.
      */
      if (n_particles_in_cell > 0)
        {
          const double n_particles_in_cell_double = static_cast<double>(n_particles_in_cell);
          const double granularity_double = static_cast<double>(granularity);
          const double ideal_n_particles_per_bucket = n_particles_in_cell_double/(Utilities::fixed_power<dim>(granularity_double));

          /*
          The "buckets_ideal" table contains doubles because the ideal
          number of particles per bucket will not always be an integer.
          */
          Table<dim,double> buckets_ideal;
          TableIndices<dim> bucket_sizes;
          for (unsigned int i=0; i<dim; ++i)
            bucket_sizes[i] = granularity;

          buckets_ideal.reinit(bucket_sizes);
          const double bucket_width = 1.0/granularity_double;
          buckets_ideal.fill(ideal_n_particles_per_bucket);

          const Table<dim,unsigned int> buckets_actual
            = sort_particles_into_buckets(particle_range, bucket_width);

          /*
          In the worst case, all particles are in one bucket.
          (granularity^dim)-1 is equal to the number of buckets
          in the table minus 1 (the bucket all the particles are in).
          */
          const double worst_case_empty_buckets = static_cast<double>(Utilities::fixed_power<dim>(granularity_double))-1.0;
          const double worst_case_error_squared =
            ((n_particles_in_cell_double-ideal_n_particles_per_bucket)*(n_particles_in_cell_double-ideal_n_particles_per_bucket))+
            ((ideal_n_particles_per_bucket*ideal_n_particles_per_bucket)*(worst_case_empty_buckets));

          double actual_error_squared = 0;
          for (unsigned int x=0; x<granularity; ++x)
            {
              for (unsigned int y=0; y<granularity; ++y)
                {
                  TableIndices<dim> entry_index;
                  entry_index[0] = x;
                  entry_index[1] = y;
                  // do another loop if in 3d
                  if (dim == 3)
                    {
                      for (unsigned int z=0; z<granularity; ++z)
                        {
                          entry_index[2] = z;
                          const double value_ideal = buckets_ideal(entry_index);
                          const double value_actual = static_cast<double>(buckets_actual(entry_index));
                          actual_error_squared += (value_ideal - value_actual)*(value_ideal - value_actual);
                        }
                    }
                  else
                    {
                      const double value_ideal = buckets_ideal(entry_index);
                      const double value_actual = static_cast<double>(buckets_actual(entry_index));
                      actual_error_squared += (value_ideal - value_actual)*(value_ideal - value_actual);
                    }
                }
            }

          return actual_error_squared/worst_case_error_squared;
        }
      else
        {
          // If there are no particles in the cell, there is no point trying to redistribute particles in the cell.
          return 0;
        }
    }



    template <int dim>
    Table<dim,unsigned int>
    ParticleHistogram<dim>::sort_particles_into_buckets(
      const typename Particles::ParticleHandler<dim>::particle_iterator_range &particle_range,
      const double bucket_width) const
    {
      TableIndices<dim> bucket_sizes;
      for (unsigned int i=0; i<dim; ++i)
        bucket_sizes[i] = granularity;
      Table<dim,unsigned int> buckets;
      buckets.reinit(bucket_sizes);
      for (const auto &particle: particle_range)
        {
          const double particle_x = particle.get_reference_location()[0];
          const double particle_y = particle.get_reference_location()[1];

          const double x_ratio = (particle_x) / (bucket_width);
          const double y_ratio = (particle_y) / (bucket_width);

          unsigned int x_index = static_cast<unsigned int>(std::floor(x_ratio));
          unsigned int y_index = static_cast<unsigned int>(std::floor(y_ratio));
          /*
          If a particle is exactly on the boundary of two cells its
          reference location will equal 1, and if this is the case,
          the "x/y/z_index" will be outside of the range of the table without
          these checks. The table has a number of entries equal to "granularity" in each dimension,
          and the table is indexed at 0, so if the "x/y/z_indez" equals "granularity" it
          will be out of range.
          */
          if (x_index == granularity)
            x_index = granularity-1;
          if (y_index == granularity)
            y_index = granularity-1;

          TableIndices<dim> entry_index;
          entry_index[0] = x_index;
          entry_index[1] = y_index;
          if (dim == 3)
            {
              const double particle_z = particle.get_reference_location()[2];
              const double z_ratio = (particle_z) / (bucket_width);
              unsigned int z_index = static_cast<unsigned int>(std::floor(z_ratio));
              if (z_index == granularity)
                z_index = granularity-1;
              entry_index[2] = z_index;
            }

          ++buckets(entry_index);
        }
      return buckets;
    }
  }
}

// explicit instantiation
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class ParticleHistogram<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
