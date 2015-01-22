/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#include <aspect/particle/generator.h>

#include <boost/random.hpp>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      // Generate random uniform distribution of particles over entire simulation domain
      template <int dim, class T>
      class RandomUniformGenerator : public Interface<dim, T>
      {
        public:
          /**
           * Constructor.
           *
           * @param[in] The MPI communicator for synchronizing particle generation.
           */
          RandomUniformGenerator() {}

          /**
           * Generate a uniformly randomly distributed set of particles in the current triangulation.
           */
          // TODO: fix the particle system so it works even with processors assigned 0 cells
          virtual
          void
          generate_particles(Particle::World<dim, T> &world,
                             const double total_num_particles)
          {
            double      total_volume, local_volume, subdomain_fraction, start_fraction, end_fraction;

            // Calculate the number of particles in this domain as a fraction of total volume
            total_volume = local_volume = 0;
            for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
                 it=world.get_triangulation()->begin_active();
                 it!=world.get_triangulation()->end(); ++it)
              {
                double cell_volume = it->measure();
                AssertThrow (cell_volume != 0, ExcMessage ("Found cell with zero volume."));

                if (it->is_locally_owned())
                  local_volume += cell_volume;
              }
            // Sum the local volumes over all nodes
            MPI_Allreduce(&local_volume, &total_volume, 1, MPI_DOUBLE, MPI_SUM, world.mpi_comm());

            // Assign this subdomain the appropriate fraction
            subdomain_fraction = local_volume/total_volume;

            // Sum the subdomain fractions so we don't miss particles from rounding and to create unique IDs
            MPI_Scan(&subdomain_fraction, &end_fraction, 1, MPI_DOUBLE, MPI_SUM, world.mpi_comm());
            start_fraction = end_fraction-subdomain_fraction;

            // Calculate start and end IDs so there are no gaps
            // TODO: this can create gaps for certain processor counts because of
            // floating point imprecision, figure out how to fix it
            const unsigned int  start_id = static_cast<unsigned int>(std::ceil(start_fraction*total_num_particles));
            const unsigned int  end_id   = static_cast<unsigned int>(fmin(std::ceil(end_fraction*total_num_particles), total_num_particles));
            const unsigned int  subdomain_particles = end_id - start_id;

            uniform_random_particles_in_subdomain(world, subdomain_particles, start_id);
          };

          /**
           * Generate a set of particles uniformly randomly distributed within the
           * specified triangulation. This is done using "roulette wheel" style
           * selection weighted by cell volume. We do cell-by-cell assignment of
           * particles because the decomposition of the mesh may result in a highly
           * non-rectangular local mesh which makes uniform particle distribution difficult.
           *
           * @param [in] world The particle world the particles will exist in
           * @param [in] num_particles The number of particles to generate in this subdomain
           * @param [in] start_id The starting ID to assign to generated particles
           */
          void uniform_random_particles_in_subdomain (Particle::World<dim, T> &world,
                                                      const unsigned int num_particles,
                                                      const unsigned int start_id)
          {
            unsigned int          i, d, v, num_tries, cur_id;
            double                total_volume, roulette_spin;
            std::map<double, LevelInd>        roulette_wheel;
            const unsigned int n_vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
            Point<dim>            pt, max_bounds, min_bounds;
            LevelInd              select_cell;

            // Create the roulette wheel based on volumes of local cells
            total_volume = 0;
            for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
                 it=world.get_triangulation()->begin_active(); it!=world.get_triangulation()->end(); ++it)
              {
                if (it->is_locally_owned())
                  {
                    // Assign an index to each active cell for selection purposes
                    total_volume += it->measure();
                    // Save the cell index and level for later access
                    roulette_wheel.insert(std::make_pair(total_volume, std::make_pair(it->level(), it->index())));
                  }
              }

            // Pick cells and assign particles at random points inside them
            cur_id = start_id;
            for (i=0; i<num_particles; ++i)
              {
                // Select a cell based on relative volume
                roulette_spin = total_volume*uniform_distribution_01(random_number_generator);
                select_cell = roulette_wheel.lower_bound(roulette_spin)->second;

                const typename parallel::distributed::Triangulation<dim>::active_cell_iterator
                it (world.get_triangulation(), select_cell.first, select_cell.second);

                // Get the bounds of the cell defined by the vertices
                for (d=0; d<dim; ++d)
                  {
                    min_bounds[d] = INFINITY;
                    max_bounds[d] = -INFINITY;
                  }
                for (v=0; v<n_vertices_per_cell; ++v)
                  {
                    pt = it->vertex(v);
                    for (d=0; d<dim; ++d)
                      {
                        min_bounds[d] = fmin(pt[d], min_bounds[d]);
                        max_bounds[d] = fmax(pt[d], max_bounds[d]);
                      }
                  }

                // Generate random points in these bounds until one is within the cell
                num_tries = 0;
                while (num_tries < 100)
                  {
                    for (d=0; d<dim; ++d)
                      {
                        pt[d] = uniform_distribution_01(random_number_generator) *
                                (max_bounds[d]-min_bounds[d]) + min_bounds[d];
                      }
                    try
                      {
                        if (it->point_inside(pt)) break;
                      }
                    catch (...)
                      {
                        // Debugging output, remove when Q4 mapping 3D sphere problem is resolved
                        //std::cerr << "Pt and cell " << pt << " " << select_cell.first << " " << select_cell.second << std::endl;
                        //for (int z=0;z<8;++z) std::cerr << "V" << z <<": " << it->vertex(z) << ", ";
                        //std::cerr << std::endl;
                        //***** MPI_Abort(communicator, 1);
                      }
                    num_tries++;
                  }
                AssertThrow (num_tries < 100, ExcMessage ("Couldn't generate particle (unusual cell shape?)."));

                // Add the generated particle to the set
                T new_particle(pt, cur_id);
                world.add_particle(new_particle, select_cell);

                cur_id++;
              }
          }

        private:
          /**
           * Random number generator and an object that describes a
           * uniform distribution on the interval [0,1]. These
           * will be used to generate particle locations at random.
           */
          boost::mt19937            random_number_generator;
          boost::uniform_01<double> uniform_distribution_01;
      };


      template <int dim, class T>
      Interface<dim,T> *
      create_generator_object (const std::string &generator_type)
      {
        if (generator_type == "random_uniform")
          return new RandomUniformGenerator<dim,T>();
        else
          Assert (false, ExcNotImplemented());

        return 0;
      }


      std::string
      generator_object_names ()
      {
        return ("random_uniform");
      }


      // explicit instantiations
      template
      Interface<2,Particle::BaseParticle<2> > *
      create_generator_object (const std::string &generator_type);
      template
      Interface<3,Particle::BaseParticle<3> > *
      create_generator_object (const std::string &generator_type);
    }
  }
}
