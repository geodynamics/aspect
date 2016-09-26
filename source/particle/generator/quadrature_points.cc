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
        types::particle_index particle_index = 0;

        const QGauss<dim> quadrature_formula(this->get_parameters().stokes_velocity_degree + 1);

        FEValues<dim> fe_values(this->get_mapping(),
                                this->get_fe(),
                                quadrature_formula,
                                update_values |
                                update_quadrature_points |
                                update_JxW_values);

        typename Triangulation<dim>::active_cell_iterator
        cell = this->get_triangulation().begin_active(),
        endc = this->get_triangulation().end();
        for (; cell != endc; cell++)
          {
            fe_values.reinit(cell);
            std::vector<Point<dim>> quadrature_points = fe_values.get_quadrature_points();
            for (typename std::vector<Point<dim>>::const_iterator q_points_itr = quadrature_points.begin();
                 q_points_itr != quadrature_points.end(); q_points_itr++)
              {
                Point<dim> particle_position_real = (*q_points_itr);

                try
                  {
                    Point<dim> particle_position_unit = this->get_mapping().transform_real_to_unit_cell(cell, (*q_points_itr));
                    if (GeometryInfo<dim>::is_inside_unit_cell(particle_position_unit))
                      {
                        const Particle<dim> particle(particle_position_real, particle_position_unit, particle_index);
                        const types::LevelInd cell_index(cell->level(), cell->index());
                        particles.insert(std::make_pair(cell_index, particle));
                        particle_index++;
                      }
                  }
                catch (typename Mapping<dim>::ExcTransformationFailed &)
                  {
                    AssertThrow (true,
                                 ExcMessage ("Couldn't generate particle (unusual cell shape?). "));
                  }
              }
          }
      }


      template <int dim>
      void
      QuadraturePoints<dim>::declare_parameters (ParameterHandler &prm)
      {}


      template <int dim>
      void
      QuadraturePoints<dim>::parse_parameters (ParameterHandler &prm)
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
                                         "Generate particles at the quadrature points of each active"
                                         "cell of the triangulation mesh.")
    }
  }
}
