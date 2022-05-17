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
      namespace internal
      {
        bool string_to_bool(const std::string &s)
        {
          return (s == "true" || s == "yes");
        }
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

        unsigned int positions_index = 0;
        const auto &mapping = this->get_mapping();
        // The unit cell of deal.II is [0,1]^dim. The limiter needs a 'unit' cell of [-.5,.5]^dim.
        const double unit_offset = 0.5;
        std::vector<double> property_minimums(n_particle_properties, std::numeric_limits<double>::max());
        std::vector<double> property_maximums(n_particle_properties, std::numeric_limits<double>::lowest());
        for (typename ParticleHandler<dim>::particle_iterator particle = particle_range.begin();
             particle != particle_range.end(); ++particle, ++positions_index)
          {
            const auto &particle_property_value = particle->get_properties();
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              {
                if (selected_properties[property_index] == true)
                  {
                    b[property_index][positions_index] = particle_property_value[property_index];
                    if (use_linear_least_squares_limiter[property_index] == true)
                      {
                        property_minimums[property_index] = std::min(b[property_index][positions_index], property_minimums[property_index]);
                        property_maximums[property_index] = std::max(b[property_index][positions_index], property_maximums[property_index]);
                      }
                  }
              }
            Point<dim> relative_particle_position = particle->get_reference_location();

            // A is accessed by A[column][row] here since we will need to append
            // columns into the qr matrix
            A[0][positions_index] = 1;
            for (unsigned int i = 1; i < n_matrix_columns; ++i)
              {
                relative_particle_position[i - 1] -= unit_offset;
                A[i][positions_index] = relative_particle_position[i - 1];
              }
          }

        ImplicitQR <Vector<double>> qr;
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
            if (selected_properties[property_index] == true)
              {
                qr.multiply_with_QT(QTb[property_index], b[property_index]);
                qr.solve(c[property_index], QTb[property_index]);
              }
          }
        const double half_h = .5;
        for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
          {
            if (selected_properties[property_index] == true)
              {
                if (use_linear_least_squares_limiter[property_index] == true)
                  {
                    c[property_index][0] = std::max(c[property_index][0], property_minimums[property_index]);
                    c[property_index][0] = std::min(c[property_index][0], property_maximums[property_index]);

                    const double max_total_slope = std::min(c[property_index][0] - property_minimums[property_index],
                                                            property_maximums[property_index] - c[property_index][0])
                                                   / half_h;
                    double current_total_slope = 0.0;
                    for (unsigned int i = 1; i < n_matrix_columns; ++i)
                      {
                        current_total_slope += std::abs(c[property_index][i]);
                      }

                    if (current_total_slope > max_total_slope && current_total_slope > std::numeric_limits<double>::min())
                      {
                        double slope_change_ratio = max_total_slope/current_total_slope;
                        for (unsigned int i = 1; i < n_matrix_columns; ++i)
                          c[property_index][i] *= slope_change_ratio;
                      }
                  }
                std::size_t positions_index = 0;
                for (typename std::vector<Point<dim>>::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++positions_index)
                  {
                    Point<dim> relative_support_point_location = mapping.transform_real_to_unit_cell(found_cell, *itr);
                    double interpolated_value = c[property_index][0];
                    for (unsigned int i = 1; i < n_matrix_columns; ++i)
                      {
                        relative_support_point_location[i - 1] -= unit_offset;
                        interpolated_value += c[property_index][i] * relative_support_point_location[i - 1];
                      }
                    if (use_linear_least_squares_limiter[property_index] == true)
                      {
                        // If the limiter is working correctly, we should not be significantly more than machine error outside of the range of the property bounds
                        Assert(interpolated_value >= property_minimums[property_index] - std::max(std::abs(property_minimums[property_index]), std::abs(property_maximums[property_index])) * 10. * std::numeric_limits<double>::epsilon(), ExcInternalError());
                        Assert(interpolated_value <= property_maximums[property_index] + std::max(std::abs(property_minimums[property_index]), std::abs(property_maximums[property_index])) * 10. * std::numeric_limits<double>::epsilon(), ExcInternalError());
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
                                  "Limit the interpolation of particle properties "
                                  "onto the cell so the value of each property is no "
                                  "smaller than its minimum and no larger than its "
                                  "maximum on the particles in each cell. If more "
                                  "than one value is specified, they will be treated "
                                  "as a list.");

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
                    use_linear_least_squares_limiter = ComponentMask(n_property_components, internal::string_to_bool(split[0]));
                  }
                else if (split.size() == n_property_components - n_internal_components)
                  {
                    std::vector<bool> parsed;
                    for (const auto &component: split)
                      parsed.push_back(internal::string_to_bool(component));
                    use_linear_least_squares_limiter = ComponentMask(parsed);
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
                                            "Uses linear least squares to obtain the slopes and center of a 2D or "
                                            "3D plane from the particle positions and a particular property value "
                                            "on those particles. "
                                            "Interpolate this property onto a vector of points. If the limiter is "
                                            "enabled then it will ensure the interpolated properties do not exceed the "
                                            "range of the minimum and maximum of the values of the property on the "
                                            "particles. Note that deal.II must be configured with BLAS and LAPACK to "
                                            "support this operation.")
    }
  }
}
