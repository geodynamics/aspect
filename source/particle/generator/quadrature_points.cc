/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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

#include <aspect/particle/generator/quadrature_points.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      QuadraturePoints<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        const QGauss<dim> quadrature_formula(this->get_parameters().stokes_velocity_degree + 1);

        types::particle_index n_particles_to_generate = quadrature_formula.size() * this->get_triangulation().n_locally_owned_active_cells();
        types::particle_index prefix_sum = 0;
        types::particle_index particle_index = 0;

        FEValues<dim> fe_values(this->get_mapping(),
                                this->get_fe(),
                                quadrature_formula,
                                update_quadrature_points);


        MPI_Scan(&n_particles_to_generate, &prefix_sum, 1, ASPECT_PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());

        particle_index = prefix_sum - n_particles_to_generate;

        typename parallel::distributed::Triangulation<dim>::active_cell_iterator
        cell = this->get_triangulation().begin_active(), endc = this->get_triangulation().end();

        for (; cell != endc; cell++)
          {
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);
                for (unsigned int i = 0; i < quadrature_formula.size(); i++)
                  {
                    const Particle<dim> particle(fe_values.get_quadrature_points()[i], quadrature_formula.get_points()[i], particle_index);
                    const types::LevelInd cell_index(cell->level(), cell->index());
                    particles.insert(std::make_pair(cell_index, particle));
                    ++particle_index;
                  }
              }
          }
      }


      template <int dim>
      void
      QuadraturePoints<dim>::declare_parameters (ParameterHandler &)
      {}


      template <int dim>
      void
      QuadraturePoints<dim>::parse_parameters (ParameterHandler &)
      {}
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(QuadraturePoints,
                                         "quadrature points",
                                         "Generates particles at the quadrature points of each active cell of "
                                         "the triangulation. Here, Gauss quadrature of degree (velocity\\_degree + 1), "
                                         "is used similarly to the assembly of Stokes matrix.")
    }
  }
}
