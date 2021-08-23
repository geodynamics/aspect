/*
  Copyright (C) 2017 - 2021 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/bilinear_least_squares.h>
#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/qr.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {

      bool string_to_bool(const std::string &s)
      {
        return (s == "true");
      }

      ComponentMask string_to_bool(const std::vector<std::string> &list)
      {
        std::vector<bool> results(list.size(), false);
        for (std::size_t i = 0; i < list.size(); i++)
          {
            results[i] = string_to_bool(list[i]);
          }
        return ComponentMask(results);
      }


      template <int dim>
      std::vector<std::vector<double>>
                                    BilinearLeastSquares<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                                                    const std::vector<Point<dim>> &positions,
                                                                                    const ComponentMask &selected_properties,
                                                                                    const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const unsigned int n_particle_properties = particle_handler.n_properties_per_particle();

        const unsigned int property_index = selected_properties.first_selected_component(selected_properties.size());

        AssertThrow(property_index != numbers::invalid_unsigned_int,
                    ExcMessage("Internal error: the particle property interpolator was "
                               "called without a specified component to interpolate."));

        const Point<dim> approximated_cell_midpoint = std::accumulate (positions.begin(), positions.end(), Point<dim>())
                                                      / static_cast<double> (positions.size());

        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Assert(positions.size() > 0,
                   ExcMessage("The particle property interpolator was not given any "
                              "positions to evaluate the particle cell_properties at."));


            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const typename ParticleHandler<dim>::particle_iterator_range particle_range =
          particle_handler.particles_in_cell(found_cell);


        std::vector<std::vector<double>> cell_properties(positions.size(),
                                                         std::vector<double>(n_particle_properties,
                                                                             numbers::signaling_nan<double>()));

        const unsigned int n_particles = std::distance(particle_range.begin(), particle_range.end());
        const unsigned int n_matrix_columns = (dim == 2) ? 3 : 4;

        // If there are too few particles, we can not perform a least squares interpolation
        // fall back to a simpler method instead.
        if (n_particles < n_matrix_columns)
          return fallback_interpolator.properties_at_points(particle_handler,
                                                            positions,
                                                            selected_properties,
                                                            found_cell);

        // Noticed that the size of matrix A is n_particles x n_matrix_columns
        // which usually is not a square matrix. Therefore, we find the
        // least squares solution of Ac=r by solving the reduced QR factorization
        // Ac = QRc = b -> Q^TQRc = Rc =Q^Tb
        // A is a std::vector of Vectors(which are it's columns) so that we create what the ImplicitQR
        // class needs.
        std::vector<Vector<double>> A(n_matrix_columns, Vector<double>(n_particles));
        std::vector<Vector<double>> b(n_particle_properties, Vector<double>(n_particles));

        AssertThrow(n_particles >= n_matrix_columns,
                    ExcMessage("At least one cell contained no particles. The 'bilinear'"
                               "interpolation scheme does not support this case. "));
        unsigned int positions_index = 0;
        const auto &mapping = this->get_mapping();
        const double unit_offset = 0.5; // The unit cell of deal.II is [0,1]^dim. The limiter needs a 'unit' cell of [-.5,.5]^dim.
        std::vector<std::pair<double, double>> property_bounds(n_particle_properties, std::pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min())); // {min, max}
        for (typename ParticleHandler<dim>::particle_iterator particle = particle_range.begin();
             particle != particle_range.end(); ++particle, ++positions_index)
          {
            const auto &particle_property_value = particle->get_properties();
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              {
                if (selected_properties[property_index])
                  {
                    b[property_index][positions_index] = particle_property_value[property_index];
                    if (use_linear_least_squares_limiter[property_index])
                      {
                        property_bounds[property_index].first = std::min(b[property_index][positions_index], property_bounds[property_index].first);
                        property_bounds[property_index].second = std::max(b[property_index][positions_index], property_bounds[property_index].second);
                      }
                  }
              }
            Point<dim> relative_particle_position = mapping.transform_real_to_unit_cell(cell, particle->get_location());
            relative_particle_position[0] -= unit_offset;
            relative_particle_position[1] -= unit_offset;
            if (dim == 3)
              relative_particle_position[2] -= unit_offset;
            // A is accessed by A[column][row] here since we will need to append
            // columns into the qr matrix
            A[0][positions_index] = 1;
            A[1][positions_index] = relative_particle_position[0];
            A[2][positions_index] = relative_particle_position[1];
            if (dim == 3)
              A[3][positions_index] = relative_particle_position[2];
          }

        ImplicitQR<Vector<double>>qr;
        for (const auto &column : A)
          qr.append_column(column);
        // If A is rank deficent then qr.append_column will not append
        // the first column that can be written as a linear combination of
        // other columns. We check that n_matrix_columns linearly independent
        // columns were added, or we rely on the fallback_interpolator

        if (qr.size() != n_matrix_columns)
          return fallback_interpolator.properties_at_points(particle_handler,
                                                            positions,
                                                            selected_properties,
                                                            found_cell);


        std::vector<Vector<double>> QTb(n_particle_properties, Vector<double>(n_matrix_columns));
        std::vector<Vector<double>> c(n_particle_properties, Vector<double>(n_matrix_columns));
        for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
          {
            if (selected_properties[property_index])
              {
                qr.multiply_with_QT(QTb[property_index], b[property_index]);
                qr.solve(c[property_index], QTb[property_index]);
              }
          }

        const double half_h = .5;
        const double threshold = 1e-12;
        for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
          {
            if (selected_properties[property_index])
              {
                if (use_linear_least_squares_limiter[property_index])
                  {
                    // Beginning of MJG and EGP's linear least squares limiter
                    c[property_index][0] = std::max(c[property_index][0], property_bounds[property_index].first);
                    c[property_index][0] = std::min(c[property_index][0], property_bounds[property_index].second);
                    if (c[property_index][0] > property_bounds[property_index].second - threshold ||
                        c[property_index][0] < property_bounds[property_index].first + threshold)
                      {
                        c[property_index][1] = 0;
                        c[property_index][2] = 0;
                        if (dim == 3)
                          {
                            c[property_index][3] = 0;
                          }
                      }
                    else
                      {
                        double upper_bound = std::min(c[property_index][0] - property_bounds[property_index].first,
                                                      property_bounds[property_index].second - c[property_index][0]);
                        double lower_bound = std::abs(c[property_index][1]) + std::abs(c[property_index][2]);
                        if (dim == 3)
                          {
                            lower_bound += std::abs(c[property_index][3]);
                          }
                        lower_bound *= half_h;
                        if (lower_bound > upper_bound)
                          {
                            c[property_index][1] = std::copysign(1, c[property_index][1]) * std::min(std::abs(c[property_index][1]), upper_bound / half_h);
                            c[property_index][2] = std::copysign(1, c[property_index][2]) * std::min(std::abs(c[property_index][2]), upper_bound / half_h);
                            if (dim == 3)
                              {
                                c[property_index][3] = std::copysign(1, c[property_index][3]) * std::min(std::abs(c[property_index][3]), upper_bound / half_h);
                              }

                            lower_bound = std::abs(c[property_index][1]) + std::abs(c[property_index][2]);
                            if (dim == 3)
                              {
                                lower_bound += std::abs(c[property_index][3]);
                              }
                            lower_bound *= half_h;
                            if (lower_bound > upper_bound)
                              {
                                if (dim == 2)
                                  {
                                    double change_in_slope = (lower_bound-upper_bound)/(half_h * 2);
                                    c[property_index][1] = std::copysign(1, c[property_index][1]) * (std::abs(c[property_index][1]) - change_in_slope);
                                    c[property_index][2] = std::copysign(1, c[property_index][2]) * (std::abs(c[property_index][2]) - change_in_slope);
                                  }
                                else
                                  {
                                    // dim == 3
                                    double change_in_slope = (lower_bound-upper_bound)/(3*half_h);
                                    double c_1 = std::copysign(1, c[property_index][1]) * (std::abs(c[property_index][1]) - change_in_slope);
                                    double c_2 = std::copysign(1, c[property_index][2]) * (std::abs(c[property_index][2]) - change_in_slope);
                                    double c_3 = std::copysign(1, c[property_index][3]) * (std::abs(c[property_index][3]) - change_in_slope);
                                    // Stage 3 Checks and corrections
                                    if (c[property_index][1] * c_1 <= 0)
                                      {
                                        c[property_index][2] = std::copysign(1, c[property_index][2]) * (std::abs(c_2) - std::abs(c_1/2));
                                        c[property_index][3] = std::copysign(1, c[property_index][3]) * (std::abs(c_3) - std::abs(c_1/2));
                                        c[property_index][1] = 0;
                                      }
                                    else if (c[property_index][2] * c_2 <= 0)
                                      {
                                        c[property_index][1] = std::copysign(1, c[property_index][1]) * (std::abs(c_1) - std::abs(c_2/2));
                                        c[property_index][3] = std::copysign(1, c[property_index][3]) * (std::abs(c_3) - std::abs(c_2/2));
                                        c[property_index][2] = 0;

                                      }
                                    else if (c[property_index][3] * c_3 <= 0)
                                      {
                                        c[property_index][1] = std::copysign(1, c[property_index][1]) * (std::abs(c_1) - std::abs(c_3/2));
                                        c[property_index][2] = std::copysign(1, c[property_index][2]) * (std::abs(c_2) - std::abs(c_3/2));
                                        c[property_index][3] = 0;
                                      }
                                    else
                                      {
                                        c[property_index][1] = c_1;
                                        c[property_index][2] = c_2;
                                        c[property_index][3] = c_3;
                                      }
                                  }
                              }
                          }
                      }
                  }
                std::size_t positions_index = 0;
                for (typename std::vector<Point<dim>>::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++positions_index)
                  {
                    Point<dim> relative_support_point_location = mapping.transform_real_to_unit_cell(found_cell, *itr);
                    relative_support_point_location[0] -= unit_offset;
                    relative_support_point_location[1] -= unit_offset;
                    if (dim == 3)
                      {
                        relative_support_point_location[2] -= unit_offset;
                      }
                    double interpolated_value = c[property_index][0] +
                                                c[property_index][1] * relative_support_point_location[0] +
                                                c[property_index][2] * relative_support_point_location[1];
                    if (dim == 3)
                      {
                        interpolated_value += c[property_index][3] * relative_support_point_location[2];
                      }
                    if (use_linear_least_squares_limiter[property_index])
                      {
                        interpolated_value = std::min(interpolated_value, property_bounds[property_index].second);
                        interpolated_value = std::max(interpolated_value, property_bounds[property_index].first);
                      }
                    cell_properties[positions_index][property_index] = interpolated_value;
                  }

              }
          }
        return cell_properties;
      }



      template <int dim>
      void
      BilinearLeastSquares<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Bilinear least squares");
              {
                prm.declare_entry("Use linear least squares limiter", "false",
                                  Patterns::List(Patterns::Bool()),
                                  "Limit the interpolation of all particle properties "
                                  "onto the cell so the value of each property is no "
                                  "smaller than its minimum and no larger than its "
                                  "maximum on the particles in each cell. Currently "
                                  "doesn't work on spherical grids.");

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
      BilinearLeastSquares<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Bilinear least squares");
              {
                const Postprocess::Particles<dim> &particle_postprocessor =
                  this->get_postprocess_manager().template get_matching_postprocessor<const Postprocess::Particles<dim>>();
                const auto &particle_property_information = particle_postprocessor.get_particle_world().get_property_manager().get_data_info();
                const unsigned int n_property_components = particle_property_information.n_components();
                const unsigned int n_internal_components = particle_property_information.get_components_by_field_name("internal: integrator properties");
                std::vector<std::string> split = Utilities::split_string_list(prm.get("Use linear least squares limiter"));
                if (split.size() == 1)
                  {
                    use_linear_least_squares_limiter = ComponentMask(n_property_components - n_internal_components, string_to_bool(split[0]));
                  }
                else if (split.size() == n_property_components - n_internal_components)
                  {
                    use_linear_least_squares_limiter = string_to_bool(split);
                  }
                else
                  {
                    AssertThrow(false, ExcMessage("The size of 'Use linear least squares limiter' should either be 1 or the number of particle properties"));
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
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(BilinearLeastSquares,
                                            "bilinear least squares",
                                            "Uses linear least squares to obtain the slopes and center of a 2 or "
                                            "3D plane from the particle positions and a particular property value "
                                            "on those particles. "
                                            "Interpolate this property onto the support points or to initiate the "
                                            "property value on a new particle. If the limiter is enabled then it "
                                            "will ensure the interpolated properties do not exceed the "
                                            "range of the minimum and maximum of the values of the property on the particles.")
    }
  }
}
