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
       * Return the cell-wise averaged properties of all tracers of the cell containing the
       * given positions.
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
        AssertThrow(n_particles != 0,
                    ExcMessage("At least one cell contained no particles. The 'constant "
                               "average' interpolation scheme does not support this case. "));
        const unsigned int n_properties = particles.begin()->second.get_properties().size();

        const unsigned int n_coefficients = 6;

        // Always matching against the 0th compositional field.
//          const typename aspect::Simulator::AdvectionField adv_field (aspect::Simulator::AdvectionField::composition(0));
        const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.compositional_fields+1);

        FEValues<dim> fe_values(this->get_mapping(),
                                this->get_fe(),
                                quadrature_formula,
                                update_values |
                                update_quadrature_points |
                                update_JxW_values);

        fe_values.reinit(found_cell);

        std::vector<std::vector<double> > properties(positions.size());
        const std::vector<Point<dim> > quadrature_points = fe_values.get_quadrature_points();

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

        index = 0;
        for (unsigned int i = 0; i < n_properties; ++i)
          {
            std::vector<double> properties_tmp(quadrature_formula.size());
            properties_tmp.clear();

            dealii::FullMatrix<double> r(6,1);
            r = 0;
            for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
                 particle != particle_range.second; ++particle, ++index)
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
            for (typename std::vector<Point<dim> >::const_iterator itr = quadrature_points.begin(); itr != quadrature_points.end(); ++itr, ++index_positions)
              {
                Point<dim> quadrature_point = *itr;
                double interpolated_value = c(0,0) + c(1,0)*(quadrature_point[0]) + c(2,0)*(quadrature_point[1]) + c(3,0)*(quadrature_point[0] * quadrature_point[1]) +  c(4,0)*(quadrature_point[0] * quadrature_point[0]) + c(5,0)*(quadrature_point[1] * quadrature_point[1]);
                properties_tmp.push_back(interpolated_value);
              }

            double local_cell_average = 0;
            double local_cell_area = 0;
            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
              {
                local_cell_area += fe_values.JxW(q);
                local_cell_average += properties_tmp[q] * fe_values.JxW(q);
              }
            local_cell_average /= local_cell_area;

            double offset = 0;
            if (local_cell_average < global_min[i])
              offset = global_min[i] - local_cell_average;
            else if (local_cell_average > global_max[i])
              offset = global_max[i] - local_cell_average;
            else
              offset = 0;

            index_positions = 0;
            for (typename std::vector<Point<dim> >::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++index_positions)
              {
                Point<dim> support_point = *itr;
                double interpolated_value = c(0,0) + c(1,0)*(support_point[0]) + c(2,0)*(support_point[1]) + c(3,0)*(support_point[0] * support_point[1]) +  c(4,0)*(support_point[0] * support_point[0]) + c(5,0)*(support_point[1] * support_point[1]) + offset;
                properties[index_positions].push_back(interpolated_value);
              }
            properties_tmp.clear();
          }
        return properties;
      }

      template <int dim>
      void
      BiquadraticLeastSquares<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Biquadratic");
              {
                prm.declare_entry ("Global maximum",
                                   boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                                   Patterns::List(Patterns::Double ()),
                                   "The maximum global property values that will be used in the bound preserving "
                                   "limiter for the discontinuous solutions during interpolation. "
                                   "The number of the input 'Global maximum' values seperated by ',' has to be "
                                   "the same as the number of particle properties");
                prm.declare_entry ("Global minimum",
                                   boost::lexical_cast<std::string>(-std::numeric_limits<double>::max()),
                                   Patterns::List(Patterns::Double ()),
                                   "The minimum global property values that will be used in the bound preserving "
                                   "limiter for the discontinuous solutions during interpolation. "
                                   "The number of the input 'Global minimum' values seperated by ',' has to be "
                                   "the same as the number of the particle properties");
              }
              prm.leave_subsection ();
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }

      template <int dim>
      void
      BiquadraticLeastSquares<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Biquadratic");
              {
                global_max = Utilities::string_to_double
                             (Utilities::split_string_list(prm.get ("Global maximum")));
                global_min = Utilities::string_to_double
                             (Utilities::split_string_list(prm.get ("Global minimum")));
              }
              prm.leave_subsection ();
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
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
                                            "Return the average of all tracer properties in the given cell.")
    }
  }
}
