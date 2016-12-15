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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/interpolator/biquadratic_least_squares.h>
#include <aspect/simulator.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the cell-wise evaluated properties of the biquadratic least squares function
       * at the positions.
       */
      template <int dim>
      std::vector<std::vector<double> >
      BiquadraticLeastSquares<dim>::properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                                         const std::vector<Point<dim> > &positions,
                                                         const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Point<dim> approximated_cell_midpoint = positions[0];
            if (positions.size() > 1)
              {
                Tensor<1,dim> direction_to_center;
                for (unsigned int i = 1; i<positions.size()-1; ++i)
                  direction_to_center += positions[i] - positions[0];
                direction_to_center /= positions.size() - 1;
                approximated_cell_midpoint += direction_to_center;
              }

            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const types::LevelInd cell_index = std::make_pair<unsigned int, unsigned int> (found_cell->level(),found_cell->index());
        const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator,
              typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator> particle_range = particles.equal_range(cell_index);

        const unsigned int n_particles = std::distance(particle_range.first,particle_range.second);
        const unsigned int n_properties = particles.begin()->second.get_properties().size();

        std::vector<std::vector<double> > properties(positions.size());

        AssertThrow(n_particles != 0,
                    ExcMessage("At least one cell contained no particles. The 'constant "
                               "average' interpolation scheme does not support this case. "));

        const unsigned int n_coefficients = 6;
        dealii::FullMatrix<double> A(n_particles,n_coefficients);
        A = 0;

        unsigned int index = 0;
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
             particle != particle_range.second; ++particle, ++index)
          {
            const Point<dim> location = particle->second.get_location();
            A(index,0) = 1;
            A(index,1) = location[0];
            A(index,2) = location[1];
            A(index,3) = location[0] * location[1];
            A(index,4) = location[0] * location[0];
            A(index,5) = location[1] * location[1];
          }

        dealii::FullMatrix<double> B(n_coefficients, n_coefficients);
        A.Tmmult(B, A, false);
        dealii::FullMatrix<double> B_inverse(B);
        B_inverse.gauss_jordan();

        for (unsigned int i = 0; i < n_properties; ++i)
          {
            dealii::FullMatrix<double> r(6,1);
            r = 0;
            for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
                 particle != particle_range.second; ++particle)
              {
                const double particle_property = particle->second.get_properties()[i];
                const Point<dim> position = particle->second.get_location();
                r(0,0) += particle_property;
                r(1,0) += particle_property * position[0];
                r(2,0) += particle_property * position[1];
                r(3,0) += particle_property * position[0] * position[1];
                r(4,0) += particle_property * position[0] * position[0];
                r(5,0) += particle_property * position[1] * position[1];
              }

            dealii::FullMatrix<double> c(6,1);
            c = 0;
            B_inverse.mmult(c, r);

            unsigned int index_positions = 0;
            for (typename std::vector<Point<dim> >::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++index_positions)
              {
                Point<dim> support_point = *itr;
                double interpolated_value = c(0,0) + c(1,0)*(support_point[0]) + c(2,0)*(support_point[1]) + c(3,0)*(support_point[0] * support_point[1]) +  c(4,0)*(support_point[0] * support_point[0]) + c(5,0)*(support_point[1] * support_point[1]);
                properties[index_positions].push_back(interpolated_value);
              }
          }
        return properties;
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(BiquadraticLeastSquares,
                                            "biquadratic",
                                            "Interpolates particle properties onto support points using a bilinear least squares in the given cell.")
    }
  }
}
