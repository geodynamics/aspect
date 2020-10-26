/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

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

#include <aspect/postprocess/gravity_point_values.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/termination_criteria/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    GravityPointValues<dim>::GravityPointValues ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      maximum_timesteps_between_outputs (std::numeric_limits<int>::max()),
      last_output_timestep (numbers::invalid_unsigned_int),
      output_file_number (numbers::invalid_unsigned_int)
    {}



    template <int dim>
    std::pair<std::string,std::string>
    GravityPointValues<dim>::execute (TableHandler &)
    {
      AssertThrow(false, ExcNotImplemented());
      return std::pair<std::string,std::string>();
    }



    template <>
    std::pair<std::string,std::string>
    GravityPointValues<3>::execute (TableHandler &statistics)
    {
      const int dim = 3;

      // Check time to see if we output gravity
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
          last_output_timestep = this->get_timestep_number();
        }
      if ((this->get_time() < last_output_time + output_interval)
          && (this->get_timestep_number() < last_output_timestep + maximum_timesteps_between_outputs)
          && (this->get_timestep_number() != 0))
        return std::pair<std::string,std::string>();

      // Up the counter of the number of the file by one, but not in
      // the very first output step. if we run postprocessors on all
      // iterations, only increase file number in the first nonlinear iteration
      const bool increase_file_number = (this->get_nonlinear_iteration() == 0) || (!this->get_parameters().run_postprocessors_on_nonlinear_iterations);
      if (output_file_number == numbers::invalid_unsigned_int)
        output_file_number = 0;
      else if (increase_file_number)
        ++output_file_number;

      const std::string file_prefix = "gravity-" + Utilities::int_to_string (output_file_number, 5);
      const std::string filename = (this->get_output_directory()
                                    + "output_gravity/"
                                    + file_prefix);

      // Get quadrature formula and increase the degree of quadrature over the velocity
      // element degree.
      const unsigned int degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree
                                  + quadrature_degree_increase;
      const QGauss<dim> quadrature_formula (degree);
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      // Get the value of the outer radius and inner radius:
      double model_outer_radius;
      double model_inner_radius;
      if (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()))
        {
          model_outer_radius = Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>>
                               (this->get_geometry_model()).outer_radius();
          model_inner_radius = Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>>
                               (this->get_geometry_model()).inner_radius();
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()))
        {
          model_outer_radius = Plugins::get_plugin_as_type<const GeometryModel::Chunk<dim>>
                               (this->get_geometry_model()).outer_radius();
          model_inner_radius = Plugins::get_plugin_as_type<const GeometryModel::Chunk<dim>>
                               (this->get_geometry_model()).inner_radius();
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()))
        {
          model_outer_radius = Plugins::get_plugin_as_type<const GeometryModel::Sphere<dim>>
                               (this->get_geometry_model()).radius();
          model_inner_radius = 0;
        }
      else
        {
          model_outer_radius = 1;
          model_inner_radius = 0;
        }

      // Get the value of the universal gravitational constant:
      const double G = aspect::constants::big_g;

      // Storing cartesian coordinates, density and JxW at local quadrature points in a vector
      // avoids to use MaterialModel and fe_values within the loops. Because the postprocessor
      // runs in parallel, the total number of local quadrature points has to be determined:
      const unsigned int n_locally_owned_cells = (this->get_triangulation().n_locally_owned_active_cells());
      const unsigned int n_quadrature_points_per_cell = quadrature_formula.size();

      // Declare the vector 'density_JxW' to store the density at quadrature points. The
      // density and the JxW are here together for simplicity in the equation (both variables
      // only appear together):
      std::vector<double> density_JxW (n_locally_owned_cells * n_quadrature_points_per_cell);
      std::vector<double> density_anomalies_JxW (n_locally_owned_cells * n_quadrature_points_per_cell);

      // Declare the vector 'position_point' to store the position of quadrature points:
      std::vector<Point<dim> > position_point (n_locally_owned_cells * n_quadrature_points_per_cell);

      // The following loop perform the storage of the position and density * JxW values
      // at local quadrature points:
      MaterialModel::MaterialModelInputs<dim> in(quadrature_formula.size(),this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature_formula.size(),this->n_compositional_fields());
      unsigned int local_cell_number = 0;
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            const std::vector<Point<dim> > &position_point_cell = fe_values.get_quadrature_points();
            in.reinit(fe_values, cell, this->introspection(), this->get_solution(), false);
            this->get_material_model().evaluate(in, out);
            for (unsigned int q = 0; q < n_quadrature_points_per_cell; ++q)
              {
                density_JxW[local_cell_number * n_quadrature_points_per_cell + q] = out.densities[q] * fe_values.JxW(q);
                density_anomalies_JxW[local_cell_number * n_quadrature_points_per_cell + q] = (out.densities[q]-reference_density) * fe_values.JxW(q);
                position_point[local_cell_number * n_quadrature_points_per_cell + q] = position_point_cell[q];
              }
            ++local_cell_number;
          }

      // Pre-Assign the coordinates of all satellites in a vector point:
      // *** First calculate the number of satellites according to the sampling scheme:
      unsigned int n_satellites;
      if (sampling_scheme == fibonacci_spiral)
        n_satellites = n_points_spiral * n_points_radius;
      else if (sampling_scheme == map)
        n_satellites = n_points_radius * n_points_longitude * n_points_latitude;
      else if (sampling_scheme == list_of_points)
        n_satellites = longitude_list.size();
      else n_satellites = 1;

      // *** Second assign the coordinates of all satellites:
      std::vector<Point<dim> > satellites_coordinate(n_satellites);
      if (sampling_scheme == map)
        {
          unsigned int p = 0;
          for (unsigned int h=0; h < n_points_radius; ++h)
            {
              for (unsigned int i=0; i < n_points_longitude; ++i)
                {
                  for (unsigned int j=0; j < n_points_latitude; ++j)
                    {
                      if (n_points_radius > 1)
                        satellites_coordinate[p][0] = minimum_radius + ((maximum_radius - minimum_radius) / (n_points_radius - 1)) * h;
                      else
                        satellites_coordinate[p][0] = minimum_radius;
                      if (n_points_longitude > 1)
                        satellites_coordinate[p][1] = (minimum_colongitude + ((maximum_colongitude - minimum_colongitude) / (n_points_longitude - 1)) * i) * numbers::PI / 180.;
                      else
                        satellites_coordinate[p][1] = minimum_colongitude * numbers::PI / 180.;
                      if (n_points_latitude > 1)
                        satellites_coordinate[p][2] = (minimum_colatitude + ((maximum_colatitude - minimum_colatitude) / (n_points_latitude - 1)) * j) * numbers::PI / 180.;
                      else
                        satellites_coordinate[p][2] = minimum_colatitude * numbers::PI / 180.;
                      ++p;
                    }
                }
            }
        }
      if (sampling_scheme == fibonacci_spiral)
        {
          const double golden_ratio = (1. + std::sqrt(5.))/2.;
          const double golden_angle = 2. * numbers::PI * (1. - 1./golden_ratio);
          unsigned int p = 0;
          for (unsigned int h=0; h < n_points_radius; ++h)
            {
              for (unsigned int s=0; s < n_points_spiral; ++s)
                {
                  if (n_points_radius > 1)
                    satellites_coordinate[p][0] = minimum_radius + ((maximum_radius - minimum_radius) / (n_points_radius - 1)) * h;
                  else
                    satellites_coordinate[p][0] = minimum_radius;
                  satellites_coordinate[p][2] = std::acos(1. - 2. * s / (n_points_spiral - 1.));
                  satellites_coordinate[p][1] = std::fmod((s*golden_angle), 2.*numbers::PI);
                  ++p;
                }
            }
        }
      if (sampling_scheme == list_of_points)
        {
          for (unsigned int p=0; p < n_satellites; ++p)
            {
              if (radius_list.size() == 1)
                satellites_coordinate[p][0] = radius_list[0];
              else
                satellites_coordinate[p][0] = radius_list[p];
              if (longitude_list[p] < 0)
                satellites_coordinate[p][1] = (360 + longitude_list[p]) * numbers::PI / 180.;
              else
                satellites_coordinate[p][1] = (longitude_list[p]) * numbers::PI / 180.;
              satellites_coordinate[p][2] = (90 - latitude_list[p]) * numbers::PI / 180. ;
            }
        }

      // open the file on rank 0 and write the headers
      std::ofstream output;
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          output.open(filename.c_str());
          AssertThrow(output,
                      ExcMessage("Unable to open file for writing: " + filename +"."));
          output << "# 1: position_satellite_r" << '\n'
                 << "# 2: position_satellite_phi" << '\n'
                 << "# 3: position_satellite_theta" << '\n'
                 << "# 4: position_satellite_x" << '\n'
                 << "# 5: position_satellite_y" << '\n'
                 << "# 6: position_satellite_z" << '\n'
                 << "# 7: gravity_x" << '\n'
                 << "# 8: gravity_y" << '\n'
                 << "# 9: gravity_z" << '\n'
                 << "# 10: gravity_norm" << '\n'
                 << "# 11: gravity_theory" << '\n'
                 << "# 12: gravity_potential" << '\n'
                 << "# 13: gravity_potential_theory" << '\n'
                 << "# 14: gravity_anomaly_x" << '\n'
                 << "# 15: gravity_anomaly_y" << '\n'
                 << "# 16: gravity_anomaly_z" << '\n'
                 << "# 17: gravity_anomaly_norm" << '\n'
                 << "# 18: gravity_gradient_xx" << '\n'
                 << "# 19: gravity_gradient_yy" << '\n'
                 << "# 20: gravity_gradient_zz" << '\n'
                 << "# 21: gravity_gradient_xy" << '\n'
                 << "# 22: gravity_gradient_xz" << '\n'
                 << "# 23: gravity_gradient_yz" << '\n'
                 << "# 24: gravity_gradient_theory_xx" << '\n'
                 << "# 25: gravity_gradient_theory_yy" << '\n'
                 << "# 26: gravity_gradient_theory_zz" << '\n'
                 << "# 27: gravity_gradient_theory_xy" << '\n'
                 << "# 28: gravity_gradient_theory_xz" << '\n'
                 << "# 29: gravity_gradient_theory_yz" << '\n'
                 << '\n';
        }

      // This is the main loop which computes gravity acceleration, potential and
      // gradients at a point located at the spherical coordinate [r, phi, theta].
      // This loop corresponds to the 3 integrals of Newton law:
      double sum_g = 0;
      double min_g = std::numeric_limits<double>::max();
      double max_g = -std::numeric_limits<double>::max();
      double sum_g_potential = 0;
      double min_g_potential = std::numeric_limits<double>::max();
      double max_g_potential = -std::numeric_limits<double>::max();
      for (unsigned int p=0; p < n_satellites; ++p)
        {

          // The spherical coordinates are shifted into cartesian to allow simplification
          // in the mathematical equation.
          std::array<double,dim> satellite_point_coordinate;
          satellite_point_coordinate[0] = satellites_coordinate[p][0];
          satellite_point_coordinate[1] = satellites_coordinate[p][1];
          satellite_point_coordinate[2] = satellites_coordinate[p][2];
          const Point<dim> position_satellite = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(satellite_point_coordinate);

          // For each point (i.e. satellite), the fourth integral goes over cells and
          // quadrature points to get the unique distance between those, to calculate
          // gravity vector components x,y,z (in tensor), potential and gradients.
          Tensor<1,dim> local_g;
          Tensor<1,dim> local_g_anomaly;
          Tensor<2,dim> local_g_gradient;
          double local_g_potential = 0;
          local_cell_number = 0;
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                for (unsigned int q = 0; q < n_quadrature_points_per_cell; ++q)
                  {
                    const double dist = (position_satellite - position_point[local_cell_number * n_quadrature_points_per_cell + q]).norm();
                    // For gravity acceleration:
                    const double KK = - G * density_JxW[local_cell_number * n_quadrature_points_per_cell + q] / std::pow(dist,3);
                    local_g += KK * (position_satellite - position_point[local_cell_number * n_quadrature_points_per_cell + q]);
                    // For gravity anomalies:
                    const double KK_anomalies = - G * density_anomalies_JxW[local_cell_number * n_quadrature_points_per_cell + q] / std::pow(dist,3);
                    local_g_anomaly += KK_anomalies * (position_satellite - position_point[local_cell_number * n_quadrature_points_per_cell + q]);
                    // For gravity potential:
                    local_g_potential -= G * density_JxW[local_cell_number * n_quadrature_points_per_cell + q] / dist;
                    // For gravity gradient:
                    const double grad_KK = G * density_JxW[local_cell_number * n_quadrature_points_per_cell + q] / std::pow(dist,5);
                    local_g_gradient[0][0] += grad_KK * (3.0
                                                         * std::pow((position_satellite[0] - position_point[local_cell_number * n_quadrature_points_per_cell + q][0]),2)
                                                         - std::pow(dist,2));
                    local_g_gradient[1][1] += grad_KK * (3.0
                                                         * std::pow((position_satellite[1] - position_point[local_cell_number * n_quadrature_points_per_cell + q][1]),2)
                                                         - std::pow(dist,2));
                    local_g_gradient[2][2] += grad_KK * (3.0
                                                         * std::pow((position_satellite[2] - position_point[local_cell_number * n_quadrature_points_per_cell + q][2]),2)
                                                         - std::pow(dist,2));
                    local_g_gradient[0][1] += grad_KK * (3.0
                                                         * (position_satellite[0] - position_point[local_cell_number * n_quadrature_points_per_cell + q][0])
                                                         * (position_satellite[1] - position_point[local_cell_number * n_quadrature_points_per_cell + q][1]));
                    local_g_gradient[0][2] += grad_KK * (3.0
                                                         * (position_satellite[0] - position_point[local_cell_number * n_quadrature_points_per_cell + q][0])
                                                         * (position_satellite[2] - position_point[local_cell_number * n_quadrature_points_per_cell + q][2]));
                    local_g_gradient[1][2] += grad_KK * (3.0
                                                         * (position_satellite[1] - position_point[local_cell_number * n_quadrature_points_per_cell + q][1])
                                                         * (position_satellite[2] - position_point[local_cell_number * n_quadrature_points_per_cell + q][2]));
                  }
                ++local_cell_number;
              }

          // Sum local gravity components over global domain:
          const Tensor<1,dim> g          = Utilities::MPI::sum (local_g, this->get_mpi_communicator());
          const Tensor<1,dim> g_anomaly  = Utilities::MPI::sum (local_g_anomaly, this->get_mpi_communicator());
          const Tensor<2,dim> g_gradient = Utilities::MPI::sum (local_g_gradient, this->get_mpi_communicator());
          const double g_potential       = Utilities::MPI::sum (local_g_potential, this->get_mpi_communicator());

          // sum gravity components for all n_satellites:
          sum_g += g.norm();
          sum_g_potential += g_potential;
          if (g.norm() > max_g)
            max_g = g.norm();
          if (g.norm() < min_g)
            min_g = g.norm();
          if (g_potential > max_g_potential)
            max_g_potential = g_potential;
          if (g_potential < min_g_potential)
            min_g_potential = g_potential;

          // analytical solution to calculate the theoretical gravity and its derivatives
          // from a uniform density model. Can only be used if concentric density profile.
          double g_theory = 0;
          double g_potential_theory = 0;
          Tensor<2,dim> g_gradient_theory;
          if (satellites_coordinate[p][0] <= model_inner_radius)
            {
              g_theory = 0;
              g_potential_theory = 2.0 * G * numbers::PI * reference_density * (std::pow(model_inner_radius,2) - std::pow(model_outer_radius,2));
            }
          else if ((satellites_coordinate[p][0] > model_inner_radius) && (satellites_coordinate[p][0] < model_outer_radius))
            {
              g_theory = G * numbers::PI * 4./3. * reference_density * (satellites_coordinate[p][0] - (std::pow(model_inner_radius,3)
                                                                        /  std::pow(satellites_coordinate[p][0],2)));
              g_potential_theory = G * numbers::PI * 4./3. * reference_density
                                   * ((std::pow(satellites_coordinate[p][0],2)/2.0) + (std::pow(model_inner_radius,3)/satellites_coordinate[p][0]))
                                   - G * numbers::PI * 2.0 * reference_density * std::pow(model_outer_radius,2);
            }
          else
            {
              g_theory = G * numbers::PI * 4./3. * reference_density * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                         / std::pow(satellites_coordinate[p][0],2);
              g_potential_theory = - G * numbers::PI * 4./3. * reference_density * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                                   /  satellites_coordinate[p][0];
              g_gradient_theory[0][0] = -G * numbers::PI * 4./3. * reference_density
                                        * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                                        * (std::pow(satellites_coordinate[p][0], 2) - 3.0 * std::pow(position_satellite[0],2))
                                        /  std::pow(satellites_coordinate[p][0],5);
              g_gradient_theory[1][1] = -G * numbers::PI * 4./3. * reference_density
                                        * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                                        * (std::pow(satellites_coordinate[p][0], 2) - 3.0 * std::pow(position_satellite[1],2))
                                        /  std::pow(satellites_coordinate[p][0],5);
              g_gradient_theory[2][2] = -G * numbers::PI * 4./3. * reference_density
                                        * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                                        * (std::pow(satellites_coordinate[p][0], 2) - 3.0 * std::pow(position_satellite[2],2))
                                        /  std::pow(satellites_coordinate[p][0],5);
              g_gradient_theory[0][1] = -G * numbers::PI * 4./3. * reference_density
                                        * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                                        * (- 3.0 * position_satellite[0] * position_satellite[1])
                                        /  std::pow(satellites_coordinate[p][0],5);
              g_gradient_theory[0][2] = -G * numbers::PI * 4./3. * reference_density
                                        * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                                        * (- 3.0 * position_satellite[0] * position_satellite[2])
                                        /  std::pow(satellites_coordinate[p][0],5);
              g_gradient_theory[1][2] = -G * numbers::PI * 4./3. * reference_density
                                        * (std::pow(model_outer_radius,3) - std::pow(model_inner_radius,3))
                                        * (- 3.0 * position_satellite[1] * position_satellite[2])
                                        /  std::pow(satellites_coordinate[p][0],5);
            }

          // write output.
          // g_gradient is here given in eotvos E (1E = 1e-9 per square seconds):
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              output << satellites_coordinate[p][0] << ' '
                     << satellites_coordinate[p][1] *180. / numbers::PI << ' '
                     << satellites_coordinate[p][2] *180. / numbers::PI << ' '
                     << position_satellite[0] << ' '
                     << position_satellite[1] << ' '
                     << position_satellite[2] << ' '
                     << std::setprecision(precision)
                     << g << ' '
                     << g.norm() << ' '
                     << g_theory << ' '
                     << g_potential << ' '
                     << g_potential_theory << ' '
                     << g_anomaly << ' '
                     << g_anomaly.norm() << ' '
                     << g_gradient[0][0] *1e9 << ' '
                     << g_gradient[1][1] *1e9 << ' '
                     << g_gradient[2][2] *1e9 << ' '
                     << g_gradient[0][1] *1e9 << ' '
                     << g_gradient[0][2] *1e9 << ' '
                     << g_gradient[1][2] *1e9 << ' '
                     << g_gradient_theory[0][0] *1e9 << ' '
                     << g_gradient_theory[1][1] *1e9 << ' '
                     << g_gradient_theory[2][2] *1e9 << ' '
                     << g_gradient_theory[0][1] *1e9 << ' '
                     << g_gradient_theory[0][2] *1e9 << ' '
                     << g_gradient_theory[1][2] *1e9 << ' '
                     << '\n';
            }
        }

      // write quantities in the statistic file
      const std::string name2("Average gravity acceleration (m/s^2)");
      statistics.add_value (name2, sum_g/n_satellites);
      statistics.set_precision (name2, precision);
      statistics.set_scientific (name2, true);

      const std::string name3("Minimum gravity acceleration (m/s^2)");
      statistics.add_value (name3, min_g);
      statistics.set_precision (name3, precision);
      statistics.set_scientific (name3, true);

      const std::string name4("Maximum gravity acceleration (m/s^2)");
      statistics.add_value (name4, max_g);
      statistics.set_precision (name4, precision);
      statistics.set_scientific (name4, true);

      const std::string name5("Average gravity potential (m^2/s^2)");
      statistics.add_value (name5, sum_g_potential/n_satellites);
      statistics.set_precision (name5, precision);
      statistics.set_scientific (name5, true);

      const std::string name6("Minimum gravity potential (m^2/s^2)");
      statistics.add_value (name6, min_g_potential);
      statistics.set_precision (name6, precision);
      statistics.set_scientific (name6, true);

      const std::string name7("Maximum gravity potential (m^2/s^2)");
      statistics.add_value (name7, max_g_potential);
      statistics.set_precision (name7, precision);
      statistics.set_scientific (name7, true);

      // up the next time we need output:
      set_last_output_time (this->get_time());
      last_output_timestep = this->get_timestep_number();
      return std::pair<std::string, std::string> ("Writing gravity output:", filename);
    }



    template <int dim>
    void
    GravityPointValues<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
          prm.declare_entry ("Sampling scheme", "map",
                             Patterns::Selection ("map|list|list of points|fibonacci spiral"),
                             "Choose the sampling scheme. By default, the map produces a "
                             "grid of equally angled points between a minimum and maximum "
                             "radius, longitude, and latitude. A list of points contains "
                             "the specific coordinates of the satellites. The fibonacci "
                             "spiral sampling scheme produces a uniformly distributed map "
                             "on the surface of sphere defined by a minimum and/or "
                             "maximum radius.");
          prm.declare_entry ("Number points fibonacci spiral", "200",
                             Patterns::Integer (0),
                             "Parameter for the fibonacci spiral sampling scheme: "
                             "This specifies the desired number of satellites per radius "
                             "layer. The default value is 200. Note that sampling "
                             "becomes more uniform with increasing number of satellites");
          prm.declare_entry ("Quadrature degree increase", "0",
                             Patterns::Integer (0),
                             "Quadrature degree increase over the velocity element "
                             "degree may be required when gravity is calculated near "
                             "the surface or inside the model. An increase in the "
                             "quadrature element adds accuracy to the gravity "
                             "solution from noise due to the model grid.");
          prm.declare_entry ("Number points radius", "1",
                             Patterns::Integer (0),
                             "Parameter for the map sampling scheme: "
                             "This specifies the number of points along "
                             "the radius (e.g. depth profile) between a minimum and "
                             "maximum radius.");
          prm.declare_entry ("Number points longitude", "1",
                             Patterns::Integer (0),
                             "Parameter for the map sampling scheme: "
                             "This specifies the number of points along "
                             "the longitude (e.g. gravity map) between a minimum and "
                             "maximum longitude.");
          prm.declare_entry ("Number points latitude", "1",
                             Patterns::Integer (0),
                             "Parameter for the map sampling scheme: "
                             "This specifies the number of points along "
                             "the latitude (e.g. gravity map) between a minimum and "
                             "maximum latitude.");
          prm.declare_entry ("Minimum radius", "0.",
                             Patterns::Double (0.0),
                             "Parameter for the map sampling scheme: "
                             "Minimum radius may be defined in or outside the model. "
                             "Prescribe a minimum radius for a sampling coverage at a "
                             "specific height.");
          prm.declare_entry ("Maximum radius", "0.",
                             Patterns::Double (0.0),
                             "Parameter for the map sampling scheme: "
                             "Maximum radius can be defined in or outside the model.");
          prm.declare_entry ("Minimum longitude", "-180.",
                             Patterns::Double (-180.0, 180.0),
                             "Parameter for the uniform distribution sampling scheme: "
                             "Gravity may be calculated for a sets of points along "
                             "the longitude between a minimum and maximum longitude.");
          prm.declare_entry ("Minimum latitude", "-90.",
                             Patterns::Double (-90.0, 90.0),
                             "Parameter for the uniform distribution sampling scheme: "
                             "Gravity may be calculated for a sets of points along "
                             "the latitude between a minimum and maximum latitude.");
          prm.declare_entry ("Maximum longitude", "180.",
                             Patterns::Double (-180.0, 180.0),
                             "Parameter for the uniform distribution sampling scheme: "
                             "Gravity may be calculated for a sets of points along "
                             "the longitude between a minimum and maximum longitude.");
          prm.declare_entry ("Maximum latitude", "90",
                             Patterns::Double (-90.0, 90.0),
                             "Parameter for the uniform distribution sampling scheme: "
                             "Gravity may be calculated for a sets of points along "
                             "the latitude between a minimum and maximum latitude.");
          prm.declare_entry ("Reference density", "3300.",
                             Patterns::Double (0.0),
                             "Gravity anomalies may be computed using density "
                             "anomalies relative to a reference density.");
          prm.declare_entry ("Precision in gravity output", "12",
                             Patterns::Integer (1),
                             "Set the precision of gravity acceleration, potential "
                             "and gradients in the gravity output and statistics file.");
          prm.declare_entry ("List of radius", "",
                             Patterns::List (Patterns::Double (0.)),
                             "Parameter for the list of points sampling scheme: "
                             "List of satellite radius coordinates. Just specify one "
                             "radius if all points values have the same radius. If "
                             "not, make sure there are as many radius as longitude "
                             "and latitude");
          prm.declare_entry ("List of longitude", "",
                             Patterns::List (Patterns::Double(-180.0, 180.0)),
                             "Parameter for the list of points sampling scheme: "
                             "List of satellite longitude coordinates.");
          prm.declare_entry ("List of latitude", "",
                             Patterns::List (Patterns::Double(-90.0, 90.0)),
                             "Parameter for the list of points sampling scheme: "
                             "List of satellite latitude coordinates.");
          prm.declare_entry ("Time between gravity output", "1e8",
                             Patterns::Double(0.0),
                             "The time interval between each generation of "
                             "gravity output files. A value of 0 indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Time steps between gravity output", boost::lexical_cast<std::string>(std::numeric_limits<int>::max()),
                             Patterns::Integer(0,std::numeric_limits<int>::max()),
                             "The maximum number of time steps between each generation of "
                             "gravity output files.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    GravityPointValues<dim>::parse_parameters (ParameterHandler &prm)
    {
      std::string gravity_subdirectory = this->get_output_directory() + "output_gravity/";
      Utilities::create_directory (gravity_subdirectory,
                                   this->get_mpi_communicator(),
                                   true);

      end_time = prm.get_double ("End time");
      if (this->convert_output_to_years())
        end_time *= year_in_seconds;

      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
          if ( (prm.get ("Sampling scheme") == "uniform distribution") || (prm.get ("Sampling scheme") == "map") )
            sampling_scheme = map;
          else if ( (prm.get ("Sampling scheme") == "list of points") || (prm.get ("Sampling scheme") == "list") )
            sampling_scheme = list_of_points;
          else if (prm.get ("Sampling scheme") == "fibonacci spiral")
            sampling_scheme = fibonacci_spiral;
          else
            AssertThrow (false, ExcMessage ("Not a valid sampling scheme."));
          quadrature_degree_increase = prm.get_integer ("Quadrature degree increase");
          n_points_spiral     = prm.get_integer("Number points fibonacci spiral");
          n_points_radius     = prm.get_integer("Number points radius");
          n_points_longitude  = prm.get_integer("Number points longitude");
          n_points_latitude   = prm.get_integer("Number points latitude");
          minimum_radius      = prm.get_double ("Minimum radius");
          maximum_radius      = prm.get_double ("Maximum radius");
          minimum_colongitude = prm.get_double ("Minimum longitude") + 180.;
          maximum_colongitude = prm.get_double ("Maximum longitude") + 180.;
          minimum_colatitude  = prm.get_double ("Minimum latitude") + 90.;
          maximum_colatitude  = prm.get_double ("Maximum latitude") + 90.;
          reference_density   = prm.get_double ("Reference density");
          precision = prm.get_integer ("Precision in gravity output");
          radius_list    = Utilities::string_to_double(Utilities::split_string_list(prm.get("List of radius")));
          longitude_list = Utilities::string_to_double(Utilities::split_string_list(prm.get("List of longitude")));
          latitude_list  = Utilities::string_to_double(Utilities::split_string_list(prm.get("List of latitude")));
          AssertThrow (longitude_list.size() == latitude_list.size(),
                       ExcMessage ("Make sure you have the same number of point coordinates in the list sampling scheme."));
          AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()) ||
                       Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()) ||
                       Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()),
                       ExcMessage ("This postprocessor can only be used if the geometry is a sphere, spherical shell or spherical chunk."));
          AssertThrow (! this->get_material_model().is_compressible(),
                       ExcMessage("This postprocessor was only tested for incompressible models so far."));
          AssertThrow (dim==3,
                       ExcMessage("This postprocessor was only tested for 3D models."));
          output_interval = prm.get_double ("Time between gravity output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
          maximum_timesteps_between_outputs = prm.get_integer("Time steps between gravity output");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    template <class Archive>
    void GravityPointValues<dim>::serialize (Archive &ar, const unsigned int)
    {
      // This deals with having the correct behavior during checkpoint/restart cycles:
      ar &last_output_time
      & last_output_timestep
      & output_file_number
      ;
    }



    template <int dim>
    void
    GravityPointValues<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output time:
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }



    template <int dim>
    void
    GravityPointValues<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);
      oa << (*this);

      status_strings["GravityPointValues"] = os.str();
    }



    template <int dim>
    void
    GravityPointValues<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("GravityPointValues") != status_strings.end())
        {
          std::istringstream is (status_strings.find("GravityPointValues")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(GravityPointValues,
                                  "gravity calculation",
                                  "A postprocessor that computes gravity, gravity anomalies, gravity "
                                  "potential and gravity gradients for a set of points (e.g. satellites) "
                                  "in or above the model surface for either a user-defined range of "
                                  "latitudes, longitudes and radius or a list of point coordinates."
                                  "Spherical coordinates in the output file are radius, colatitude "
                                  "and colongitude. Gravity is here based on the density distribution "
                                  "from the material model (and non adiabatic). This means that the "
                                  "density may come directly from an ascii file. This postprocessor also "
                                  "computes theoretical gravity and its derivatives, which corresponds to "
                                  "the analytical solution of gravity in the same geometry but filled "
                                  "with a reference density. The reference density is also used to "
                                  "determine density anomalies for computing gravity anomalies. Thus "
                                  "one must carefully evaluate the meaning of the gravity anomaly output, "
                                  "because the solution may not reflect the actual gravity anomaly (due to "
                                  "differences in the assumed reference density). On way to guarantee correct "
                                  "gravity anomalies is to subtract gravity of a certain point from the average "
                                  "gravity on the map. Another way is to directly use density anomalies for this "
                                  "postprocessor."
                                  "The average- minimum- and maximum gravity acceleration and potential are "
                                  "written into the statistics file.")
  }
}
