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

#include <aspect/particle/generator/random_uniform.h>

#include <boost/random.hpp>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      RandomUniform<dim>::RandomUniform() {}

      template <int dim>
      std::multimap<types::LevelInd, Particle<dim> >
      RandomUniform<dim>::generate_particles()
      {
        // Calculate the number of particles in this domain as a fraction of total volume
        double local_volume = 0;
        for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
             it=this->get_triangulation().begin_active();
             it!=this->get_triangulation().end(); ++it)
          {
            const double cell_volume = it->measure();
            AssertThrow (cell_volume != 0, ExcMessage ("Found cell with zero volume."));

            if (it->is_locally_owned())
              local_volume += cell_volume;
          }

        // Sum the local volumes over all nodes
        const double total_volume = Utilities::MPI::sum(local_volume,this->get_mpi_communicator());

        // Assign this subdomain the appropriate fraction
        double subdomain_fraction = local_volume / total_volume;

        // Sum the subdomain fractions so we don't miss particles from rounding and to create unique IDs
        double end_fraction = 0.0;
        MPI_Scan(&subdomain_fraction, &end_fraction, 1, MPI_DOUBLE, MPI_SUM, this->get_mpi_communicator());
        const double start_fraction = end_fraction - subdomain_fraction;

        // Calculate start and end IDs so there are no gaps
        // TODO: this can create gaps for certain processor counts because of
        // floating point imprecision, figure out how to fix it
        const types::particle_index start_id = static_cast<types::particle_index>(std::ceil(start_fraction*n_tracers));
        const types::particle_index end_id   = static_cast<types::particle_index>(fmin(std::ceil(end_fraction*n_tracers), n_tracers));
        const types::particle_index subdomain_particles = end_id - start_id;

        return uniform_random_particles_in_subdomain(subdomain_particles, start_id);
      }

      template <int dim>
      std::multimap<types::LevelInd, Particle<dim> >
      RandomUniform<dim>::uniform_random_particles_in_subdomain (const types::particle_index subdomain_particles,
                                                                 const types::particle_index start_id)
      {
        std::map<double, types::LevelInd> roulette_wheel;
        std::multimap<types::LevelInd, Particle<dim> > particles;

        // Create the roulette wheel based on volumes of local cells
        double total_volume = 0;
        for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
             it=this->get_triangulation().begin_active(); it!=this->get_triangulation().end(); ++it)
          {
            if (it->is_locally_owned())
              {
                // Assign an index to each active cell for selection purposes
                total_volume += it->measure();
                // Save the cell index and level for later access
                roulette_wheel.insert(std::make_pair(total_volume, types::LevelInd(it->level(), it->index())));
              }
          }

        // Pick cells and assign particles at random points inside them
        for (types::particle_index cur_id = start_id; cur_id < start_id + subdomain_particles; ++cur_id)
          {
            // Select a cell based on relative volume
            const double roulette_spin = total_volume*uniform_distribution_01(random_number_generator);
            const types::LevelInd select_cell = roulette_wheel.lower_bound(roulette_spin)->second;

            const typename parallel::distributed::Triangulation<dim>::active_cell_iterator
            it (&(this->get_triangulation()), select_cell.first, select_cell.second);

            Point<dim> pt, max_bounds, min_bounds;
            // Get the bounds of the cell defined by the vertices
            for (unsigned int d=0; d<dim; ++d)
              {
                min_bounds[d] = INFINITY;
                max_bounds[d] = -INFINITY;
              }
            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
              {
                pt = it->vertex(v);
                for (unsigned int d=0; d<dim; ++d)
                  {
                    min_bounds[d] = fmin(pt[d], min_bounds[d]);
                    max_bounds[d] = fmax(pt[d], max_bounds[d]);
                  }
              }

            // Generate random points in these bounds until one is within the cell
            unsigned int num_tries = 0;
            while (num_tries < 100)
              {
                for (unsigned int d=0; d<dim; ++d)
                  {
                    pt[d] = uniform_distribution_01(random_number_generator) *
                            (max_bounds[d]-min_bounds[d]) + min_bounds[d];
                  }
                try
                  {
                    const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(it, pt);
                    if (GeometryInfo<dim>::is_inside_unit_cell(p_unit)) break;
                  }
                catch (typename Mapping<dim>::ExcTransformationFailed &)
                  {
                    // The point is not in this cell. Do nothing, just try again.
                  }
                num_tries++;
              }
            AssertThrow (num_tries < 100, ExcMessage ("Couldn't generate particle (unusual cell shape?)."));

            // Add the generated particle to the set
            const Particle<dim> new_particle(pt, cur_id);
            particles.insert(std::make_pair(select_cell,new_particle));
          }
        return particles;
      }


      template <int dim>
      void
      RandomUniform<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            prm.declare_entry ("Number of tracers", "1000",
                               Patterns::Double (0),
                               "Total number of tracers to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      RandomUniform<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            n_tracers    = static_cast<unsigned int>(prm.get_double ("Number of tracers"));
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(RandomUniform,
                                         "random uniform",
                                         "Generate random uniform distribution of "
                                         "particles over entire simulation domain.")
    }
  }
}
