/*
  Copyright (C) 2017 - 2024 by the authors of the ASPECT code.

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
#include <aspect/particle/manager.h>
#include <aspect/utilities.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/qr.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
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

        const typename ParticleHandler<dim>::particle_iterator_range particle_range =
          particle_handler.particles_in_cell(cell);

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
                                                            cell);

        // Noticed that the size of matrix A is n_particles x n_matrix_columns
        // which usually is not a square matrix. Therefore, we find the
        // least squares solution of Ac=r by solving the reduced QR factorization
        // Ac = QRc = b -> Q^TQRc = Rc =Q^Tb
        // A is a std::vector of Vectors(which are it's columns) so that we
        // create what the ImplicitQR class needs.
        std::vector<Vector<double>> A(n_matrix_columns, Vector<double>(n_particles));
        std::vector<Vector<double>> b(n_particle_properties, Vector<double>(n_particles));

        unsigned int particle_index = 0;
        // The unit cell of deal.II is [0,1]^dim. The limiter needs a 'unit' cell of [-.5,.5]^dim.
        const double unit_offset = 0.5;
        std::vector<double> property_minimums(n_particle_properties, std::numeric_limits<double>::max());
        std::vector<double> property_maximums(n_particle_properties, std::numeric_limits<double>::lowest());
        for (typename ParticleHandler<dim>::particle_iterator particle = particle_range.begin();
             particle != particle_range.end(); ++particle, ++particle_index)
          {
            const ArrayView<double> particle_property_value = particle->get_properties();
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              {
                if (selected_properties[property_index] == true)
                  {
                    b[property_index][particle_index] = particle_property_value[property_index];
                    if (use_linear_least_squares_limiter[property_index] == true)
                      {
                        property_minimums[property_index] = std::min(property_minimums[property_index], particle_property_value[property_index]);
                        property_maximums[property_index] = std::max(property_maximums[property_index], particle_property_value[property_index]);
                      }
                  }
              }
            Point<dim> relative_particle_position = particle->get_reference_location();

            // A is accessed by A[column][row] here since we will need to append
            // columns into the qr matrix
            A[0][particle_index] = 1;
            for (unsigned int i = 1; i < n_matrix_columns; ++i)
              {
                relative_particle_position[i - 1] -= unit_offset;
                A[i][particle_index] = relative_particle_position[i - 1];
              }
          }

        // If the limiter is enabled for at least one property then we know that we can access ghost cell
        // particles to determine the bounds of the properties on the mode (due to the assert of
        // 'Exchange ghost particles' in parse_parameters). Otherwise we do not need to access those particles
        if (use_linear_least_squares_limiter.n_selected_components() != 0)
          {
            std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> active_neighbors;
            GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim>>(cell, active_neighbors);
            for (const auto &active_neighbor : active_neighbors)
              {
                if (active_neighbor->is_artificial())
                  continue;
                const std::vector<double> neighbor_cell_average = fallback_interpolator.properties_at_points(particle_handler, positions, selected_properties, active_neighbor)[0];
                for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
                  {
                    if (selected_properties[property_index] == true && use_linear_least_squares_limiter[property_index] == true)
                      {
                        property_minimums[property_index] = std::min(property_minimums[property_index], neighbor_cell_average[property_index]);
                        property_maximums[property_index] = std::max(property_maximums[property_index], neighbor_cell_average[property_index]);
                      }
                  }
              }
            if (cell->at_boundary())
              {
                const std::vector<double> cell_average_values = fallback_interpolator.properties_at_points(particle_handler, {positions[0]}, selected_properties, cell)[0];
                for (unsigned int face_id = 0; face_id < cell->reference_cell().n_faces(); ++face_id)
                  {
                    if (cell->at_boundary(face_id))
                      {
                        const unsigned int opposing_face_id = GeometryInfo<dim>::opposite_face[face_id];
                        const auto &opposing_cell = cell->neighbor(opposing_face_id);
                        if (opposing_cell.state() == IteratorState::IteratorStates::valid && opposing_cell->is_active() && opposing_cell->is_artificial() == false)
                          {
                            const auto neighbor_cell_average = fallback_interpolator.properties_at_points(particle_handler, {positions[0]}, selected_properties, opposing_cell)[0];
                            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
                              {
                                if (selected_properties[property_index] == true && use_boundary_extrapolation[property_index] == true)
                                  {
                                    Assert(cell->reference_cell().is_hyper_cube() == true, ExcNotImplemented());
                                    const double expected_boundary_value = 1.5 * cell_average_values[property_index] - 0.5 * neighbor_cell_average[property_index];
                                    property_minimums[property_index] = std::min(property_minimums[property_index], expected_boundary_value);
                                    property_maximums[property_index] = std::max(property_maximums[property_index], expected_boundary_value);
                                  }
                              }
                          }
                      }
                  }
              }
          }

        ImplicitQR<Vector<double>> qr;
        for (const auto &column : A)
          qr.append_column(column);
        // If A is rank deficent then qr.append_column will not append
        // the first column that can be written as a linear combination of
        // other columns. We check that all columns were added or we
        // rely on the fallback interpolator
        if (qr.size() != n_matrix_columns)
          return fallback_interpolator.properties_at_points(particle_handler,
                                                            positions,
                                                            selected_properties,
                                                            cell);

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
                    Point<dim> relative_support_point_location = this->get_mapping().transform_real_to_unit_cell(cell, *itr);
                    double interpolated_value = c[property_index][0];
                    for (unsigned int i = 1; i < n_matrix_columns; ++i)
                      {
                        relative_support_point_location[i - 1] -= unit_offset;
                        interpolated_value += c[property_index][i] * relative_support_point_location[i - 1];
                      }
                    if (use_linear_least_squares_limiter[property_index] == true)
                      {
                        // Assert that the limiter was reasonably effective. We can not expect perfect accuracy
                        // due to inaccuracies e.g. in the inversion of the mapping.
                        const double tolerance = std::sqrt(std::numeric_limits<double>::epsilon())
                                                 * std::max(std::abs(property_minimums[property_index]),
                                                            std::abs(property_maximums[property_index]));
                        (void) tolerance;
                        Assert(interpolated_value >= property_minimums[property_index] - tolerance,
                               ExcMessage("The particle interpolation limiter did not succeed. Interpolated value: " + std::to_string(interpolated_value)
                                          + " is smaller than the minimum particle property value: " + std::to_string(property_minimums[property_index]) + "."));
                        Assert(interpolated_value <= property_maximums[property_index] + tolerance,
                               ExcMessage("The particle interpolation limiter did not succeed. Interpolated value: " + std::to_string(interpolated_value)
                                          + " is larger than the maximum particle property value: " + std::to_string(property_maximums[property_index]) + "."));

                        interpolated_value = std::min(interpolated_value, property_maximums[property_index]);
                        interpolated_value = std::max(interpolated_value, property_minimums[property_index]);
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
        prm.enter_subsection("Interpolator");
        {
          prm.enter_subsection("Bilinear least squares");
          {
            prm.declare_entry("Use linear least squares limiter", "true",
                              Patterns::List(Patterns::Bool()),
                              "Limit the interpolation of particle properties onto the cell, so that "
                              "the value of each property is no smaller than its minimum and no "
                              "larger than its maximum on the particles of each cell, and the "
                              "average of neighboring cells. If more than one value is given, "
                              "it will be treated as a list with one component per particle property.");
            prm.declare_entry("Use boundary extrapolation", "false",
                              Patterns::List(Patterns::Bool()),
                              "Extends the range used by 'Use linear least squares limiter' "
                              "by linearly interpolating values at cell boundaries from neighboring "
                              "cells. If more than one value is given, it will be treated as a list "
                              "with one component per particle property. Enabling 'Use boundary "
                              "extrapolation' requires enabling 'Use linear least squares "
                              "limiter'.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      BilinearLeastSquares<dim>::parse_parameters (ParameterHandler &prm)
      {
        fallback_interpolator.parse_parameters(prm);
        prm.enter_subsection("Interpolator");
        {
          prm.enter_subsection("Bilinear least squares");
          {
            const auto &particle_property_information = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager().get_data_info();
            const unsigned int n_property_components = particle_property_information.n_components();
            const unsigned int n_internal_components = particle_property_information.get_components_by_field_name("internal: integrator properties");

            std::vector<std::string> linear_least_squares_limiter_split = Utilities::split_string_list(prm.get("Use linear least squares limiter"));
            std::vector<bool> linear_least_squares_limiter_parsed;
            if (linear_least_squares_limiter_split.size() == 1)
              {
                linear_least_squares_limiter_parsed = std::vector<bool>(n_property_components - n_internal_components, Utilities::string_to_bool(linear_least_squares_limiter_split[0]));
              }
            else if (linear_least_squares_limiter_split.size() == n_property_components - n_internal_components)
              {
                for (const auto &component: linear_least_squares_limiter_split)
                  linear_least_squares_limiter_parsed.push_back(Utilities::string_to_bool(component));
              }
            else
              {
                AssertThrow(false, ExcMessage("The size of 'Use linear least squares limiter' should either be 1 or the number of particle properties"));
              }
            for (unsigned int i = 0; i < n_internal_components; ++i)
              linear_least_squares_limiter_parsed.push_back(false);
            use_linear_least_squares_limiter = ComponentMask(linear_least_squares_limiter_parsed);



            const std::vector<std::string> boundary_extrapolation_split = Utilities::split_string_list(prm.get("Use boundary extrapolation"));
            std::vector<bool> boundary_extrapolation_parsed;
            if (boundary_extrapolation_split.size() == 1)
              {
                boundary_extrapolation_parsed = std::vector<bool>(n_property_components - n_internal_components, Utilities::string_to_bool(boundary_extrapolation_split[0]));
              }
            else if (boundary_extrapolation_split.size() == n_property_components - n_internal_components)
              {
                for (const auto &component: boundary_extrapolation_split)
                  boundary_extrapolation_parsed.push_back(Utilities::string_to_bool(component));
              }
            else
              {
                AssertThrow(false, ExcMessage("The size of 'Use boundary extrapolation' should either be 1 or the number of particle properties"));
              }
            for (unsigned int i = 0; i < n_internal_components; ++i)
              boundary_extrapolation_parsed.push_back(false);
            use_boundary_extrapolation = ComponentMask(boundary_extrapolation_parsed);
            for (unsigned int property_index = 0; property_index < n_property_components - n_internal_components; ++property_index)
              {
                AssertThrow(use_linear_least_squares_limiter[property_index] || !use_boundary_extrapolation[property_index],
                            ExcMessage("'Use boundary extrapolation' must be set with 'Use linear least squares limiter' to be valid."));
              }
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
                                            "Uses linear least squares to obtain the slopes and center of a 2d or "
                                            "3d plane from the particle positions and a particular property value "
                                            "on those particles. "
                                            "Interpolate this property onto a vector of points. If the limiter is "
                                            "enabled then it will ensure the interpolated properties do not exceed the "
                                            "range of the minimum and maximum of the values of the property on the "
                                            "particles. Note that deal.II must be configured with BLAS and LAPACK to "
                                            "support this operation.")
    }
  }
}
