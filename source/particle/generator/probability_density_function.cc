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

#include <aspect/particle/generator/probability_density_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/geometry_info.h>

#include <boost/random.hpp>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      ProbabilityDensityFunction<dim>::ProbabilityDensityFunction()
        :
        random_number_generator(5432),
        n_tracers()
      {}

      template <int dim>
      std::multimap<types::LevelInd, Particle<dim> >
      ProbabilityDensityFunction<dim>::generate_particles()
      {
        // Get the local accumulated probabilities for every cell
        const std::vector<double> accumulated_cell_weights = local_accumulated_cell_weights();

        // Sum the local integrals over all nodes
        double local_function_integral = accumulated_cell_weights.back();
        const double global_function_integral = Utilities::MPI::sum (local_function_integral,
                                                                     this->get_mpi_communicator());

        AssertThrow(global_function_integral > std::numeric_limits<double>::min(),
                    ExcMessage("The integral of the user prescribed probability density function over the domain equals zero, "
                               "ASPECT has no way to determine the cell of generated tracers. Please ensure that the provided "
                               "function is positive in at least a part of the domain, also check the syntax of the function."));

        // Determine the starting weight of this process, which is the sum of
        // the weights of all processes with a lower rank
        double start_weight = 0.0;
        MPI_Scan(&local_function_integral, &start_weight, 1, MPI_DOUBLE, MPI_SUM, this->get_mpi_communicator());
        start_weight -= local_function_integral;

        // Build the map between integrated probability density and cell
        std::map<double, types::LevelInd> cell_weights;
        typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active(),
                                                       endc = this->get_dof_handler().end();
        for (unsigned int i = 0; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              // adjust the cell_weights according to the local starting weight
              cell_weights.insert(std::make_pair(accumulated_cell_weights[i] + start_weight, types::LevelInd(cell->level(),cell->index())));
              ++i;
            }

        // Calculate start and end IDs so there are no gaps
        const types::particle_index start_id = llround(n_tracers * start_weight / global_function_integral) + 1;

        return generate_particles_in_subdomain(cell_weights,global_function_integral,start_weight,n_tracers, 1.1 * start_id);
      }

      template <int dim>
      std::vector<double>
      ProbabilityDensityFunction<dim>::local_accumulated_cell_weights () const
      {
        //evaluate function at all cell midpoints, sort cells according to weight
        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_quadrature_points = quadrature_formula.size();

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_quadrature_points |
                                 update_JxW_values);

        std::vector<double> accumulated_cell_weights;

        // compute the integral weight by quadrature
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              const std::vector<Point<dim> > position = fe_values.get_quadrature_points();

              // Then add the weight of the current cell. Weights are always
              // interpreted positively, even if the function evaluates to a
              // negative number.
              double next_cell_weight = 0.0;
              if (accumulated_cell_weights.size() > 0)
                next_cell_weight = accumulated_cell_weights.back();

              for (unsigned int q = 0; q < n_quadrature_points; ++q)
                next_cell_weight += std::fabs(function.value(position[q]) * fe_values.JxW(q));

              // Start from the weight of the previous cell
              accumulated_cell_weights.push_back(next_cell_weight);
            }
        return accumulated_cell_weights;
      }

      template <int dim>
      std::multimap<types::LevelInd, Particle<dim> >
      ProbabilityDensityFunction<dim>::generate_particles_in_subdomain (const std::map<double,types::LevelInd> &cells,
                                                                        const double global_weight,
                                                                        const double start_weight,
                                                                        const types::particle_index num_particles,
                                                                        const types::particle_index start_id)
      {
        std::multimap<types::LevelInd, Particle<dim> > particles;

        // Pick cells and assign particles at random points inside them
        for (types::particle_index cur_id = start_id; cur_id < start_id + num_particles; ++cur_id)
          {
            // Draw the random number that determines the cell of the tracer
            const double random_weight =  global_weight * uniform_distribution_01(random_number_generator);

            // If the selected cell is not owned by this process, continue with next tracer
            const std::map<double,types::LevelInd>::const_iterator select_cell_iterator = cells.lower_bound(random_weight);
            if ((random_weight < start_weight) || (select_cell_iterator == cells.end()))
              continue;

            const types::LevelInd select_cell = select_cell_iterator->second;

            const typename parallel::distributed::Triangulation<dim>::active_cell_iterator
            it (&(this->get_triangulation()), select_cell.first, select_cell.second);

            Point<dim> max_bounds, min_bounds;
            // Get the bounds of the cell defined by the vertices
            for (unsigned int d=0; d<dim; ++d)
              {
                min_bounds[d] = INFINITY;
                max_bounds[d] = -INFINITY;
              }
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
              {
                const Point<dim> vertex_position = it->vertex(v);
                for (unsigned int d=0; d<dim; ++d)
                  {
                    min_bounds[d] = fmin(vertex_position[d], min_bounds[d]);
                    max_bounds[d] = fmax(vertex_position[d], max_bounds[d]);
                  }
              }

            // Generate random points in these bounds until one is within the cell
            unsigned int num_tries = 0;
            Point<dim> particle_position;
            while (num_tries < 100)
              {
                for (unsigned int d=0; d<dim; ++d)
                  {
                    particle_position[d] = uniform_distribution_01(random_number_generator) *
                                           (max_bounds[d]-min_bounds[d]) + min_bounds[d];
                  }
                try
                  {
                    const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(it, particle_position);
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
            const Particle<dim> new_particle(particle_position, cur_id);
            particles.insert(std::make_pair(select_cell, new_particle));
          }

        return particles;
      }



      template <int dim>
      void
      ProbabilityDensityFunction<dim>::declare_parameters (ParameterHandler &prm)
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

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Random function");
              {
                Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
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
      ProbabilityDensityFunction<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            n_tracers    = static_cast<types::particle_index>(prm.get_double ("Number of tracers"));

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Random function");
              {
                try
                  {
                    function.parse_parameters (prm);
                  }
                catch (...)
                  {
                    std::cerr << "ERROR: FunctionParser failed to parse\n"
                              << "\t'Initial conditions.Function'\n"
                              << "with expression\n"
                              << "\t'" << prm.get("Function expression") << "'";
                    throw;
                  }
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(ProbabilityDensityFunction,
                                         "probability density function",
                                         "Generate random distribution of "
                                         "particles over entire simulation domain. "
                                         "The particle density is prescribed in the "
                                         "form of a user-prescribed function. The "
                                         "format of this function follows the syntax "
                                         "understood by the muparser library, see "
                                         "Section~\\ref{sec:muparser-format}. The "
                                         "return value of the function is always "
                                         "interpreted as a positive probability "
                                         "density but it can be set to zero in"
                                         "parts of the domain.")
    }
  }
}
