/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/quadratic_least_squares.h>
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
      double QuadraticLeastSquares<dim>::evaluate_interpolation_function(const Vector<double> &coefficients, const Point<dim> &position) const
      {
        if (dim == 2)
          {
            return coefficients[0] +
                   coefficients[1] * position[0] +
                   coefficients[2] * position[1] +
                   coefficients[3] * position[0] * position[1] +
                   coefficients[4] * position[0] * position[0] +
                   coefficients[5] * position[1] * position[1];
          }
        else
          {
            return coefficients[0] +
                   coefficients[1] * position[0] +
                   coefficients[2] * position[1] +
                   coefficients[3] * position[2] +
                   coefficients[4] * position[0] * position[1] +
                   coefficients[5] * position[0] * position[2] +
                   coefficients[6] * position[1] * position[2] +
                   coefficients[7] * position[0] * position[0] +
                   coefficients[8] * position[1] * position[1] +
                   coefficients[9] * position[2] * position[2];
          }
      }


      template <int dim>
      std::pair<double, double> QuadraticLeastSquares<dim>::get_interpolation_bounds(const Vector<double> &coefficients) const
      {
        double interpolation_min = std::numeric_limits<double>::max();
        double interpolation_max = std::numeric_limits<double>::lowest();
        for (const auto &critical_point : get_critical_points(coefficients))
          {
            bool critical_point_in_cell = true;
            for (unsigned int d = 0; d < dim; ++d)
              {
                if (critical_point[d] < -0.5 || critical_point[d] > 0.5)
                  critical_point_in_cell = false;
              }
            if (critical_point_in_cell)
              {
                const double value_at_critical_point = evaluate_interpolation_function(coefficients, critical_point);
                interpolation_min = std::min(interpolation_min, value_at_critical_point);
                interpolation_max = std::max(interpolation_max, value_at_critical_point);
              }
          }
        return {interpolation_min, interpolation_max};
      }


      template <int dim>
      std::vector<Point<dim>> QuadraticLeastSquares<dim>::get_critical_points(const Vector<double> &coefficients) const
      {
        std::vector<Point<dim>> critical_points;
        const double epsilon = 10. * coefficients.linfty_norm() * std::numeric_limits<double>::epsilon();
        if (dim == 2)
          {
            // reserve the maximum number of critical points
            // in 2d: one inside, 4 edges, and 4 corners
            critical_points.reserve(1 + 4 + 4);
            // If finding the critical point of the function (or along a cell edge) would
            // require division by 0, or the solve of a singular matrix, then there is not
            // a unique critical point. There cannot be two critical points to this function,
            // so there must be infinitely many points with the same value, so it should be
            // caught by one of the other checks of the value of the interpolation

            // compute the location of the global critical point
            Tensor<2, dim, double> critical_point_A;
            critical_point_A[0][0] = 2 * coefficients[4];
            critical_point_A[0][1] = coefficients[3];
            critical_point_A[1][0] = coefficients[3];
            critical_point_A[1][1] = 2 * coefficients[5];
            Tensor<1, dim, double> critical_point_b;
            critical_point_b[0] = -coefficients[1];
            critical_point_b[1] = -coefficients[2];
            if (std::abs(determinant(critical_point_A)) > epsilon)
              {
                critical_points.emplace_back(invert(critical_point_A) * critical_point_b);
              }

            // Compute the critical value for each of the edges. This is necessary even if we found the
            // critical point inside the unit cell, because the value at the edges can be a minimum, while
            // the critical point inside the cell is a maximum, or vice-versa. Additionally the critical
            // point could be a saddle point, in which case we would still need to find a minimum and maximum over the cell.
            if (std::abs(coefficients[5]) > epsilon)
              {
                critical_points.emplace_back(-0.5, -(2 * coefficients[2] - coefficients[3])/(4 * coefficients[5]));
                critical_points.emplace_back( 0.5, -(2 * coefficients[2] + coefficients[3])/(4 * coefficients[5]));
              }
            if (std::abs(coefficients[4]) > epsilon)
              {
                critical_points.emplace_back(-(2 * coefficients[1] - coefficients[3])/(4 * coefficients[4]), -0.5);
                critical_points.emplace_back(-(2 * coefficients[1] + coefficients[3])/(4 * coefficients[4]),  0.5);
              }

            // Compute the critical value for each of the corners. This is necessary even if critical points
            // have already been found in previous steps, as the global critical point could be a minimum,
            // and the edge critical points could also be minimums.
            for (double x = -0.5; x <= 0.5; ++x)
              {
                for (double y = -0.5; y <= 0.5; ++y)
                  {
                    critical_points.emplace_back(x,y);
                  }
              }
          }
        else if (dim == 3)
          {
            // reserve the maximum number of critical points
            // in 3d: one inside, 6 faces, 12 edges, and 8 corners
            critical_points.reserve(1 + 6 + 12 + 8);
            // If finding the critical point of the function (or along a cell edge) would
            // require division by 0, or the solve of a singular matrix, then there is not
            // a unique critical point. There cannot be two critical points to this function,
            // so there must be infinitely many points with the same value, so it should be
            // caught by one of the other checks of the value of the interpolation

            // Compute the location of the global critical point
            {
              Tensor<2, dim, double> critical_point_A;
              critical_point_A[0][0] = 2 * coefficients[7];
              critical_point_A[0][1] = coefficients[4];
              critical_point_A[0][2] = coefficients[5];
              critical_point_A[1][0] = coefficients[4];
              critical_point_A[1][1] = 2 * coefficients[8];
              critical_point_A[1][2] = coefficients[6];
              critical_point_A[2][0] = coefficients[5];
              critical_point_A[2][1] = coefficients[6];
              critical_point_A[2][2] = 2 * coefficients[9];
              Tensor<1, dim, double> critical_point_b;
              critical_point_b[0] = -coefficients[1];
              critical_point_b[1] = -coefficients[2];
              critical_point_b[2] = -coefficients[3];
              if (std::abs(determinant(critical_point_A)) > epsilon)
                {
                  critical_points.emplace_back(invert(critical_point_A) * critical_point_b);
                }
            }

            // Compute the location of critical points along the faces of the cell.
            // This is is necessary even if we found a global critical point as it
            // could be a minimum and the faces could have a maximum or vice-versa.
            Tensor<2, 2, double> critical_point_A;
            Tensor<1, 2, double> critical_point_b;
            Tensor<1, 2, double> critical_point_X;
            // The columns of this critical_point_A correspond to Y and Z.
            critical_point_A[0][0] = 2 * coefficients[8];
            critical_point_A[0][1] = coefficients[6];
            critical_point_A[1][0] = coefficients[6];
            critical_point_A[1][1] = 2 * coefficients[9];
            if (std::abs(determinant(critical_point_A)) > epsilon)
              {
                const Tensor<2, 2, double> critical_point_A_inv = invert(critical_point_A);
                double x = -0.5;
                critical_point_b[0] = -(coefficients[2] + coefficients[4] * x);
                critical_point_b[1] = -(coefficients[3] + coefficients[5] * x);
                critical_point_X = critical_point_A_inv * critical_point_b;
                critical_points.emplace_back(x, critical_point_X[0], critical_point_X[1]);
                x = 0.5;
                critical_point_b[0] = -(coefficients[2] + coefficients[4] * x);
                critical_point_b[1] = -(coefficients[3] + coefficients[5] * x);
                critical_point_X = critical_point_A_inv * critical_point_b;
                critical_points.emplace_back(x, critical_point_X[0], critical_point_X[1]);
              }
            // The columns of this critical_point_A correspond to X and Z.
            critical_point_A[0][0] = 2 * coefficients[7];
            critical_point_A[0][1] = coefficients[5];
            critical_point_A[1][0] = coefficients[5];
            critical_point_A[1][1] = 2 * coefficients[9];
            if (std::abs(determinant(critical_point_A)) > epsilon)
              {
                const Tensor<2, 2, double> critical_point_A_inv = invert(critical_point_A);
                double y = -0.5;
                critical_point_b[0] = -(coefficients[1] + coefficients[4] * y);
                critical_point_b[1] = -(coefficients[3] + coefficients[6] * y);
                critical_point_X = critical_point_A_inv * critical_point_b;
                critical_points.emplace_back(critical_point_X[0], y, critical_point_X[1]);
                y = 0.5;
                critical_point_b[0] = -(coefficients[1] + coefficients[4] * y);
                critical_point_b[1] = -(coefficients[3] + coefficients[6] * y);
                critical_point_X = critical_point_A_inv * critical_point_b;
                critical_points.emplace_back(critical_point_X[0], y, critical_point_X[1]);
              }
            // The columns of this critical_point_A correspond to X and Y.
            critical_point_A[0][0] = 2 * coefficients[7];
            critical_point_A[0][1] = coefficients[4];
            critical_point_A[1][0] = coefficients[4];
            critical_point_A[1][1] = 2 * coefficients[8];
            if (std::abs(determinant(critical_point_A)) > epsilon)
              {
                const Tensor<2, 2, double> critical_point_A_inv = invert(critical_point_A);
                double z = -0.5;
                critical_point_b[0] = -(coefficients[1] + coefficients[5] * z);
                critical_point_b[1] = -(coefficients[2] + coefficients[6] * z);
                critical_point_X = critical_point_A_inv * critical_point_b;
                critical_points.emplace_back(critical_point_X[0], critical_point_X[1], z);
                z = 0.5;
                critical_point_b[0] = -(coefficients[1] + coefficients[5] * z);
                critical_point_b[1] = -(coefficients[2] + coefficients[6] * z);
                critical_point_X = critical_point_A_inv * critical_point_b;
                critical_points.emplace_back(critical_point_X[0], critical_point_X[1], z);
              }

            // Compute the location of critical points along the edges.
            // This is necessary even if critical points have been found in previous
            // steps, as the global critial point and critical points on faces could
            // all be minimums.
            if (std::abs(coefficients[9]) > epsilon)
              {
                for (double x = -0.5; x <= 0.5; ++x)
                  {
                    for (double y = -0.5; y <= 0.5; ++y)
                      {
                        critical_points.emplace_back(x,y, -(coefficients[3] + coefficients[5] * x + coefficients[6] * y)/(2 * coefficients[9]));
                      }
                  }
              }
            if (std::abs(coefficients[8]) > epsilon)
              {
                for (double x = -0.5; x <= 0.5; ++x)
                  {
                    for (double z = -0.5; z <= 0.5; ++z)
                      {
                        critical_points.emplace_back(x, -(coefficients[2] + coefficients[4] * x + coefficients[6] * z) / (2 * coefficients[8]), z);
                      }
                  }
              }
            if (std::abs(coefficients[7]) > epsilon)
              {
                for (double y = -0.5; y <= 0.5; ++y)
                  {
                    for (double z = -0.5; z <= 0.5; ++z)
                      {
                        critical_points.emplace_back(-(coefficients[1] + coefficients[4] * y + coefficients[5] * z)/(2*coefficients[7]), y, z);
                      }
                  }
              }

            // Compute the location of critical points along the corners
            // This is necessary even if critical points have been found in previous
            // steps, as the previously found critical points could all be minimums
            // and the corners could hold the maximum value over the cell.
            for (double x = -0.5; x <= 0.5; ++x)
              {
                for (double y = -0.5; y <= 0.5; ++y)
                  {
                    for (double z = -0.5; z <= 0.5; ++z)
                      {
                        critical_points.emplace_back(x, y, z);
                      }
                  }
              }
          }
        return critical_points;
      }

      template <int dim>
      std::vector<std::vector<double>>
      QuadraticLeastSquares<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
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

        const unsigned int n_matrix_columns = (dim == 2) ? 6 : 10;
        if (n_particles < n_matrix_columns)
          return fallback_interpolator.properties_at_points(particle_handler,
                                                            positions,
                                                            selected_properties,
                                                            cell);
        const std::vector<double> cell_average_values = fallback_interpolator.properties_at_points(particle_handler,
        {positions[0]},
        selected_properties,
        cell)[0];


        // Notice that the size of matrix A is n_particles x n_matrix_columns
        // which usually is not a square matrix. Therefore, we find the
        // least squares solution of Ac=r by solving the reduced QR factorization
        // Ac = QRc = b -> Q^TQRc = Rc =Q^Tb
        // A is a std::vector of Vectors(which are it's columns) so that we
        // create what the ImplicitQR class needs.
        std::vector<Vector<double>> A(n_matrix_columns, Vector<double>(n_particles));
        std::vector<Vector<double>> b(n_particle_properties, Vector<double>(n_particles));

        unsigned int particle_index = 0;
        // The unit cell of deal.II is [0, 1]^dim. The limiter needs a 'unit' cell of [-0.5, 0.5]^dim
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
                    property_minimums[property_index] = std::min(property_minimums[property_index], particle_property_value[property_index]);
                    property_maximums[property_index] = std::max(property_maximums[property_index], particle_property_value[property_index]);
                  }
              }

            Point<dim> relative_particle_position = particle->get_reference_location();
            for (unsigned int d = 0; d < dim; ++d)
              relative_particle_position[d] -= unit_offset;
            // A is accessed by A[column][row] here since we need to append
            // columns into the qr matrix.

            // There is a potential that in the future we may change to
            // interpolate to $Q_2$ instead of the current $P_2$. This would
            // involve adding more terms to the interpolation for some problems.
            // We for now leave those terms out of the interpolation because they
            // would require more particles per cell and we are not yet aware of
            // a benchmark where those terms have affected the convegence.
            A[0][particle_index] = 1;
            A[1][particle_index] = relative_particle_position[0];
            A[2][particle_index] = relative_particle_position[1];
            if (dim == 2)
              {
                A[3][particle_index] = relative_particle_position[0] * relative_particle_position[1];
                A[4][particle_index] = relative_particle_position[0] * relative_particle_position[0];
                A[5][particle_index] = relative_particle_position[1] * relative_particle_position[1];
              }
            else
              {
                A[3][particle_index] = relative_particle_position[2];
                A[4][particle_index] = relative_particle_position[0] * relative_particle_position[1];
                A[5][particle_index] = relative_particle_position[0] * relative_particle_position[2];
                A[6][particle_index] = relative_particle_position[1] * relative_particle_position[2];
                A[7][particle_index] = relative_particle_position[0] * relative_particle_position[0];
                A[8][particle_index] = relative_particle_position[1] * relative_particle_position[1];
                A[9][particle_index] = relative_particle_position[2] * relative_particle_position[2];
              }
          }

        // If the limiter is enabled for at least one property then we know that we can access ghost cell
        // particles to determine the bounds of the properties on the model (due to the assert of
        // 'Exchange ghost particles' in parse_parameters). Otherwise we do not need to access those particles
        if (use_quadratic_least_squares_limiter.n_selected_components(n_particle_properties) != 0)
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
                    if (selected_properties[property_index] == true && use_quadratic_least_squares_limiter[property_index] == true)
                      {
                        property_minimums[property_index] = std::min(property_minimums[property_index], neighbor_cell_average[property_index]);
                        property_maximums[property_index] = std::max(property_maximums[property_index], neighbor_cell_average[property_index]);
                      }
                  }

              }
            if (cell->at_boundary())
              {
                for (unsigned int face_id = 0; face_id < cell->reference_cell().n_faces(); ++face_id)
                  {
                    if (cell->at_boundary(face_id))
                      {
                        const unsigned int opposing_face_id = GeometryInfo<dim>::opposite_face[face_id];
                        const auto &opposing_cell = cell->neighbor(opposing_face_id);
                        if (opposing_cell.state() == IteratorState::IteratorStates::valid && opposing_cell->is_active() && !opposing_cell->is_artificial())
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
            if (selected_properties[property_index])
              {
                qr.multiply_with_QT(QTb[property_index], b[property_index]);
                qr.solve(c[property_index], QTb[property_index]);
                if (use_quadratic_least_squares_limiter[property_index])
                  {

                    const std::pair<double, double> interpolation_bounds = get_interpolation_bounds(c[property_index]);
                    const double interpolation_min = interpolation_bounds.first;
                    const double interpolation_max = interpolation_bounds.second;
                    if ((interpolation_max - cell_average_values[property_index]) > std::numeric_limits<double>::epsilon() &&
                        (cell_average_values[property_index] - interpolation_min) > std::numeric_limits<double>::epsilon())
                      {
                        const double alpha = std::max(std::min((cell_average_values[property_index] - property_minimums[property_index])/(cell_average_values[property_index] - interpolation_min),
                                                               (property_maximums[property_index]-cell_average_values[property_index])/(interpolation_max - cell_average_values[property_index])), 0.0);
                        // If alpha > 1, then using it would make the function grow to meet the bounds.
                        if (alpha < 1.0)
                          {
                            c[property_index] *= alpha;
                            c[property_index][0] += (1-alpha) * cell_average_values[property_index];
                          }
                      }
                  }
              }
          }
        unsigned int index_positions = 0;
        for (typename std::vector<Point<dim>>::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++index_positions)
          {
            Point<dim> relative_support_point_location = this->get_mapping().transform_real_to_unit_cell(cell, *itr);
            for (unsigned int d = 0; d < dim; ++d)
              relative_support_point_location[d] -= unit_offset;
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              {
                if (selected_properties[property_index] == true)
                  {
                    double interpolated_value = evaluate_interpolation_function(c[property_index], relative_support_point_location);
                    // Overshoot and undershoot correction of interpolated particle property.
                    if (use_quadratic_least_squares_limiter[property_index])
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

                        // This chopping is done to avoid values that are just outside
                        // of the limiting bounds.
                        interpolated_value = std::min(interpolated_value, property_maximums[property_index]);
                        interpolated_value = std::max(interpolated_value, property_minimums[property_index]);
                      }
                    cell_properties[index_positions][property_index] = interpolated_value;
                  }

              }
          }
        return cell_properties;
      }



      template <int dim>
      void
      QuadraticLeastSquares<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Interpolator");
        {
          prm.enter_subsection("Quadratic least squares");
          {
            prm.declare_entry("Use quadratic least squares limiter", "true",
                              Patterns::List(Patterns::Bool()),
                              "Limit the interpolation of particle properties onto the cell, so that "
                              "the value of each property is no smaller than its minimum and no "
                              "larger than its maximum on the particles of each cell, and the "
                              "average of neighboring cells. If more than one value is given, "
                              "it will be treated as a list with one component per particle property.");
            prm.declare_entry("Use boundary extrapolation", "false",
                              Patterns::List(Patterns::Bool()),
                              "Extends the range used by 'Use quadratic least squares limiter' "
                              "by linearly interpolating values at cell boundaries from neighboring "
                              "cells. If more than one value is given, it will be treated as a list "
                              "with one component per particle property. Enabling 'Use boundary "
                              "extrapolation' requires enabling 'Use quadratic least squares "
                              "limiter'.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      QuadraticLeastSquares<dim>::parse_parameters (ParameterHandler &prm)
      {
        fallback_interpolator.parse_parameters(prm);

        prm.enter_subsection("Interpolator");
        {
          prm.enter_subsection("Quadratic least squares");
          {
            const auto &particle_property_information = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager().get_data_info();
            const unsigned int n_property_components = particle_property_information.n_components();
            const unsigned int n_internal_components = particle_property_information.get_components_by_field_name("internal: integrator properties");

            const std::vector<std::string> quadratic_least_squares_limiter_split = Utilities::split_string_list(prm.get("Use quadratic least squares limiter"));
            std::vector<bool> quadratic_least_squares_limiter_parsed;
            if (quadratic_least_squares_limiter_split.size() == 1)
              {
                quadratic_least_squares_limiter_parsed = std::vector<bool>(n_property_components - n_internal_components, Utilities::string_to_bool(quadratic_least_squares_limiter_split[0]));
              }
            else if (quadratic_least_squares_limiter_split.size() == n_property_components - n_internal_components)
              {
                for (const auto &component: quadratic_least_squares_limiter_split)
                  quadratic_least_squares_limiter_parsed.push_back(Utilities::string_to_bool(component));
              }
            else
              {
                AssertThrow(false, ExcMessage("The size of 'Use quadratic least squares limiter' should either be 1 or the number of particle properties"));
              }
            for (unsigned int i = 0; i < n_internal_components; ++i)
              quadratic_least_squares_limiter_parsed.push_back(false);
            use_quadratic_least_squares_limiter = ComponentMask(quadratic_least_squares_limiter_parsed);


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
                AssertThrow(use_quadratic_least_squares_limiter[property_index] || !use_boundary_extrapolation[property_index],
                            ExcMessage("'Use boundary extrapolation' must be set with 'Use quadratic least squares limiter' to be valid."));
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
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(QuadraticLeastSquares,
                                            "quadratic least squares",
                                            "Interpolates particle properties onto a vector of points using a "
                                            "quadratic least squares method. Note that deal.II must be configured "
                                            "with BLAS/LAPACK.")
    }
  }
}
