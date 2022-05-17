/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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
      void QuadraticLeastSquares<dim>::update_bounds(const Vector<double> &coefficients, const Point<dim> &position, double &interpolation_min, double &interpolation_max) const
      {
        for (unsigned int d = 0; d < dim; ++d)
          if (position[d] < -0.5 || position[d] > 0.5)
            return;
        const double value_at_position = evaluate_interpolation_function(coefficients, position);
        interpolation_min = std::min(interpolation_min, value_at_position);
        interpolation_max = std::max(interpolation_max, value_at_position );
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
                              "positions to evaluate the particle properties at."));


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

        const unsigned int n_matrix_columns = (dim == 2) ? 6 : 10;
        if (n_particles < n_matrix_columns)
          return fallback_interpolator.properties_at_points(particle_handler,
                                                            positions,
                                                            selected_properties,
                                                            found_cell);
        const std::vector<double> cell_average_values = fallback_interpolator.properties_at_points(particle_handler,
        {positions[0]},
        selected_properties,
        found_cell)[0];


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
            const auto &particle_property_value = particle->get_properties();
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              {
                if (selected_properties[property_index] == true)
                  {
                    b[property_index][particle_index] = particle_property_value[property_index];
                    property_minimums[property_index] = std::min(property_minimums[property_index], particle_property_value[property_index]);
                    property_minimums[property_index] = std::min(property_minimums[property_index], particle_property_value[property_index]);
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
        std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> active_neighbors;
        GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim>>(found_cell, active_neighbors);
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
        if (found_cell->at_boundary())
          {
            for (unsigned int face_id = 0; face_id < found_cell->reference_cell().n_faces(); ++face_id)
              {
                if (found_cell->at_boundary(face_id))
                  {
                    const unsigned int opposing_face_id = GeometryInfo<dim>::opposite_face[face_id];
                    const auto &opposing_cell = found_cell->neighbor(opposing_face_id);
                    if (!opposing_cell->is_locally_owned())
                      {
                        continue;
                      }
                    const auto neighbor_cell_average = fallback_interpolator.properties_at_points(particle_handler, {positions[0]}, selected_properties, opposing_cell)[0];
                    for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
                      {
                        if (selected_properties[property_index] == true && use_boundary_extrapolation[property_index] == true)
                          {
                            Assert(found_cell->reference_cell().is_hyper_cube() == true, ExcNotImplemented());
                            const double expected_boundary_value = 1.5 * cell_average_values[property_index] - 0.5 * neighbor_cell_average[property_index];
                            property_minimums[property_index] = std::min(property_minimums[property_index], expected_boundary_value);
                            property_maximums[property_index] = std::max(property_maximums[property_index], expected_boundary_value);
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
                                                            found_cell);
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
                    double interpolation_min = std::numeric_limits<double>::max();
                    double interpolation_max = std::numeric_limits<double>::lowest();
                    if (dim == 2)
                      {
                        // If finding the critical point of the function (or along a cell edge) would
                        // require division by 0, then there is not a unique critical point or it is a
                        // saddle point. There cannot be two critical points to this function, so there
                        // must be infinitely many with the same value, so it should be caught by one of
                        // the other checks of the value of the function

                        // compute the location of the critical point
                        Point<dim> critical_point;
                        Tensor<2, dim, double> critical_point_A;
                        critical_point_A[0][0] = 2 * c[property_index][4];
                        critical_point_A[0][1] = c[property_index][3];
                        critical_point_A[1][0] = c[property_index][3];
                        critical_point_A[1][1] = 2 * c[property_index][5];
                        Tensor<1, dim, double> critical_point_b;
                        critical_point_b[0] = -c[property_index][1];
                        critical_point_b[1] = -c[property_index][2];
                        if (std::abs(determinant(critical_point_A)) > std::numeric_limits<double>::epsilon())
                          {
                            critical_point = invert(critical_point_A) * critical_point_b;
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                          }

                        // Compute the critical value for each of the edges. This is necessary even if we found the
                        // critical point inside the unit cell, because the value at the edges can be a minimum, while
                        // the critical point inside the cell is a maximum, or vice-versa. Additionally the critical
                        // point could be a saddle point, in which case we would still need to find a minimum and maximum over the cell.
                        if (std::abs(c[property_index][5]) > std::numeric_limits<double>::epsilon())
                          {
                            critical_point[0] = -0.5;
                            critical_point[1] = -(2 * c[property_index][2] - c[property_index][3])/(4 * c[property_index][5]);
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);

                            critical_point[0] = 0.5;
                            critical_point[1] = -(2 * c[property_index][2] + c[property_index][3])/(4 * c[property_index][5]);
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                          }
                        if (std::abs(c[property_index][4]) > std::numeric_limits<double>::epsilon())
                          {
                            critical_point[0] = -(2 * c[property_index][1] - c[property_index][3])/(4 * c[property_index][4]);
                            critical_point[1] = -0.5;
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);

                            critical_point[0] = -(2 * c[property_index][1] + c[property_index][3])/(4 * c[property_index][4]);
                            critical_point[1] = 0.5;
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                          }
                        // Compute the critical value for each of the corners. This is neccessary even if we
                        // have found critical points in previous steps, as the
                        for (double x = -0.5; x <= 0.5; ++x)
                          {
                            for (double y = -0.5; y <= 0.5; ++y)
                              {
                                critical_point[0] = x;
                                critical_point[1] = y;
                                update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                              }
                          }
                      }
                    else
                      {
                        // If finding the critical point of the function (or along a cell edge) would
                        // require division by 0, then there is not a unique critical point. There
                        // cannot be two critical points to this function, so there must be infinitely
                        // many with the same value, so it should be caught by one of the other checks
                        // of the value of the function
                        Point<dim> critical_point;
                        {
                          Tensor<2, dim, double> critical_point_A;
                          critical_point_A[0][0] = 2*c[property_index][7];
                          critical_point_A[0][1] = c[property_index][4];
                          critical_point_A[0][2] = c[property_index][5];
                          critical_point_A[1][0] = c[property_index][4];
                          critical_point_A[1][1] = 2 * c[property_index][8];
                          critical_point_A[1][2] = c[property_index][6];
                          critical_point_A[2][0] = c[property_index][5];
                          critical_point_A[2][1] = c[property_index][6];
                          critical_point_A[2][2] = 2 * c[property_index][9];
                          Tensor<1, dim, double> critical_point_b;
                          critical_point_b[0] = -c[property_index][1];
                          critical_point_b[1] = -c[property_index][2];
                          critical_point_b[2] = -c[property_index][3];
                          if (std::abs(determinant(critical_point_A)) > std::numeric_limits<double>::epsilon())
                            {
                              critical_point = invert(critical_point_A) * critical_point_b;
                              update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                            }
                        }
                        // Faces

                        Tensor<2, 2, double> critical_point_A;
                        Tensor<1, 2, double> critical_point_b;
                        Tensor<1, 2, double> critical_point_X;
                        // The columns of this critical_point_A correspond to Y and Z.
                        critical_point_A[0][0] = 2 * c[property_index][8];
                        critical_point_A[0][1] = c[property_index][6];
                        critical_point_A[1][0] = c[property_index][6];
                        critical_point_A[1][1] = 2 * c[property_index][9];
                        if (std::abs(determinant(critical_point_A)) > std::numeric_limits<double>::epsilon())
                          {
                            critical_point[0] = -0.5;
                            critical_point_b[0] = -(c[property_index][2] + c[property_index][4] * critical_point[0]);
                            critical_point_b[1] = -(c[property_index][3] + c[property_index][5] * critical_point[0]);
                            critical_point_X = invert(critical_point_A) * critical_point_b;
                            critical_point[1] = critical_point_X[0];
                            critical_point[2] = critical_point_X[1];
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                            critical_point[0] = 0.5;
                            critical_point_b[0] = -(c[property_index][2] + c[property_index][4] * critical_point[0]);
                            critical_point_b[1] = -(c[property_index][3] + c[property_index][5] * critical_point[0]);
                            critical_point_X = invert(critical_point_A) * critical_point_b;
                            critical_point[1] = critical_point_X[0];
                            critical_point[2] = critical_point_X[1];
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                          }
                        // The columns of this critical_point_A correspond to X and Z.
                        critical_point_A[0][0] = 2 * c[property_index][7];
                        critical_point_A[0][1] = c[property_index][5];
                        critical_point_A[1][0] = c[property_index][5];
                        critical_point_A[1][1] = 2 * c[property_index][9];
                        if (std::abs(determinant(critical_point_A)) > std::numeric_limits<double>::epsilon())
                          {
                            critical_point[1] = -0.5;
                            critical_point_b[0] = -(c[property_index][1] + c[property_index][4] * critical_point[1]);
                            critical_point_b[1] = -(c[property_index][3] + c[property_index][6] * critical_point[1]);
                            critical_point_X = invert(critical_point_A) * critical_point_b;
                            critical_point[0] = critical_point_X[0];
                            critical_point[2] = critical_point_X[1];
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                            critical_point[1] = 0.5;
                            critical_point_b[0] = -(c[property_index][1] + c[property_index][4] * critical_point[1]);
                            critical_point_b[1] = -(c[property_index][3] + c[property_index][6] * critical_point[1]);
                            critical_point_X = invert(critical_point_A) * critical_point_b;
                            critical_point[0] = critical_point_X[0];
                            critical_point[2] = critical_point_X[1];
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                          }
                        // The columns of this critical_point_A correspond to X and Y.
                        critical_point_A[0][0] = 2 * c[property_index][7];
                        critical_point_A[0][1] = c[property_index][4];
                        critical_point_A[1][0] = c[property_index][4];
                        critical_point_A[1][1] = 2 * c[property_index][8];
                        if (std::abs(determinant(critical_point_A)) > std::numeric_limits<double>::epsilon())
                          {
                            critical_point[2] = -0.5;
                            critical_point_b[0] = -(c[property_index][1] + c[property_index][5] * critical_point[2]);
                            critical_point_b[1] = -(c[property_index][2] + c[property_index][6] * critical_point[2]);
                            critical_point_X = invert(critical_point_A) * critical_point_b;
                            critical_point[0] = critical_point_X[0];
                            critical_point[1] = critical_point_X[1];
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                            critical_point[2] = 0.5;
                            critical_point_b[0] = -(c[property_index][1] + c[property_index][5] * critical_point[2]);
                            critical_point_b[1] = -(c[property_index][2] + c[property_index][6] * critical_point[2]);
                            critical_point_X = invert(critical_point_A) * critical_point_b;
                            critical_point[0] = critical_point_X[0];
                            critical_point[1] = critical_point_X[1];
                            update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                          }

                        // Edges
                        if (std::abs(c[property_index][9]) > std::numeric_limits<double>::epsilon())
                          {
                            for (double x = -0.5; x <= 0.5; ++x)
                              {
                                for (double y = -0.5; y <= 0.5; ++y)
                                  {
                                    critical_point[0] = x;
                                    critical_point[1] = y;
                                    critical_point[2] = -(c[property_index][3] + c[property_index][5] * critical_point[0] + c[property_index][6] * critical_point[1])/(2 * c[property_index][9]);
                                    update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                                  }
                              }
                          }
                        if (std::abs(c[property_index][7]) > std::numeric_limits<double>::epsilon())
                          {
                            for (double y = -0.5; y <= 0.5; ++y)
                              {
                                for (double z = -0.5; z <= 0.5; ++z)
                                  {
                                    critical_point[1] = y;
                                    critical_point[2] = z;
                                    critical_point[0] = -(c[property_index][1] + c[property_index][4] * critical_point[1] + c[property_index][5] * critical_point[2])/(2*c[property_index][7]);
                                    update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                                  }
                              }
                          }
                        if (std::abs(c[property_index][8]) > std::numeric_limits<double>::epsilon())
                          {
                            for (double x = -0.5; x <= 0.5; ++x)
                              {
                                for (double z = -0.5; z <= 0.5; ++z)
                                  {
                                    critical_point[0] = x;
                                    critical_point[2] = z;
                                    critical_point[1] = -(c[property_index][1] + c[property_index][4] * critical_point[0] + c[property_index][6] * critical_point[2]) / (2 * c[property_index][8]);
                                    update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                                  }
                              }
                          }
                        // Corners
                        for (double x = -0.5; x <= 0.5; ++x)
                          {
                            for (double y = -0.5; y <= 0.5; ++y)
                              {
                                for (double z = -0.5; z <= 0.5; ++z)
                                  {
                                    critical_point[0] = x;
                                    critical_point[1] = y;
                                    critical_point[2] = z;
                                    update_bounds(c[property_index], critical_point, interpolation_min, interpolation_max);
                                  }
                              }
                          }
                      }
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
            Point<dim> relative_support_point_location = this->get_mapping().transform_real_to_unit_cell(found_cell, *itr);
            for (unsigned int d = 0; d < dim; ++d)
              relative_support_point_location[d] -= unit_offset;
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              {
                double interpolated_value = evaluate_interpolation_function(c[property_index], relative_support_point_location);
                // Overshoot and undershoot correction of interpolated particle property.
                if (use_quadratic_least_squares_limiter[property_index])
                  {
                    Assert(interpolated_value >= property_minimums[property_index] - std::max(std::abs(property_minimums[property_index]), std::abs(property_maximums[property_index])) * 10. * std::numeric_limits<double>::epsilon(), ExcInternalError());
                    Assert(interpolated_value <= property_maximums[property_index] + std::max(std::abs(property_minimums[property_index]), std::abs(property_maximums[property_index])) * 10. * std::numeric_limits<double>::epsilon(), ExcInternalError());

                    interpolated_value = std::min(interpolated_value, property_maximums[property_index]);
                    interpolated_value = std::max(interpolated_value, property_minimums[property_index]);
                  }

                cell_properties[index_positions][property_index] = interpolated_value;
              }
          }
        return cell_properties;
      }



      template <int dim>
      void
      QuadraticLeastSquares<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Quadratic least squares");
              {
                prm.declare_entry("Use quadratic least squares limiter", "true",
                                  Patterns::List(Patterns::Bool()),
                                  "Limit the interpolation of particle properties onto the cell, so that"
                                  "the value of each property is no smaller than its minimum and no "
                                  "larger than its maximum on the particles of each cell, and the "
                                  "average of neighboring cells. If more than one value is given, "
                                  "it will be treated as a list with one component per particle property.");
                prm.declare_entry("Use boundary extrapolation", "false",
                                  Patterns::List(Patterns::Bool()),
                                  "Extends the range used by 'Use quadratic limiter' by linearly "
                                  "interpolating values at cell boundaries from neighboring cells. "
                                  "If more than one value is given, it will be treated as a list with one component per particle property. "
                                  "Enabling 'Use boundary extrapolation' without also enabling "
                                  "'Use quadratic least squares limiter' does nothing.");
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
      QuadraticLeastSquares<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Quadratic least squares");
              {
                const Postprocess::Particles<dim> &particle_postprocessor =
                  this->get_postprocess_manager().template get_matching_postprocessor<const Postprocess::Particles<dim>>();
                const auto &particle_property_information = particle_postprocessor.get_particle_world().get_property_manager().get_data_info();
                const unsigned int n_property_components = particle_property_information.n_components();
                const unsigned int n_internal_components = particle_property_information.get_components_by_field_name("internal: integrator properties");

                const std::vector<std::string> quadratic_least_squares_limiter_split = Utilities::split_string_list(prm.get("Use quadratic least squares limiter"));
                if (quadratic_least_squares_limiter_split.size() == 1)
                  {
                    use_quadratic_least_squares_limiter = ComponentMask(n_property_components, internal::string_to_bool(quadratic_least_squares_limiter_split[0]));
                  }
                else if (quadratic_least_squares_limiter_split.size() == n_property_components - n_internal_components)
                  {
                    std::vector<bool> parsed;
                    for (const auto &component: quadratic_least_squares_limiter_split)
                      parsed.push_back(internal::string_to_bool(component));
                    use_quadratic_least_squares_limiter = ComponentMask(parsed);
                  }
                else
                  {
                    AssertThrow(false, ExcMessage("The size of 'Use quadratic least squares limiter' should either be 1 or the number of particle properties"));
                  }
                const std::vector<std::string> boundary_extrapolation_split = Utilities::split_string_list(prm.get("Use boundary extrapolation"));
                if (boundary_extrapolation_split.size() == 1)
                  {
                    use_boundary_extrapolation = ComponentMask(n_property_components, internal::string_to_bool(boundary_extrapolation_split[0]));
                  }
                else if (boundary_extrapolation_split.size() == n_property_components - n_internal_components)
                  {
                    std::vector<bool> parsed;
                    for (const auto &component: boundary_extrapolation_split)
                      parsed.push_back(internal::string_to_bool(component));
                    use_boundary_extrapolation = ComponentMask(parsed);
                  }
                else
                  {
                    AssertThrow(false, ExcMessage("The size of 'Use boundary extrapolation' should either be 1 or the number of particle properties"));
                  }
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
