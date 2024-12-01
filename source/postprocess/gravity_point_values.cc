/*
  Copyright (C) 2018 - 2024 by the authors of the ASPECT code.

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
#include <deal.II/base/mpi.templates.h>
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
    void
    GravityPointValues<dim>::initialize ()
    {
      // Get the value of the outer radius and inner radius:
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
      satellite_positions_spherical.resize(n_satellites);
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
                        satellite_positions_spherical[p][0] = minimum_radius + ((maximum_radius - minimum_radius) / (n_points_radius - 1)) * h;
                      else
                        satellite_positions_spherical[p][0] = minimum_radius;
                      if (n_points_longitude > 1)
                        satellite_positions_spherical[p][1] = (minimum_colongitude + ((maximum_colongitude - minimum_colongitude) /
                                                                                      (n_points_longitude - 1)) * i) * constants::degree_to_radians;
                      else
                        satellite_positions_spherical[p][1] = minimum_colongitude * constants::degree_to_radians;
                      if (n_points_latitude > 1)
                        satellite_positions_spherical[p][2] = (minimum_colatitude + ((maximum_colatitude - minimum_colatitude) /
                                                                                     (n_points_latitude - 1)) * j) * constants::degree_to_radians;
                      else
                        satellite_positions_spherical[p][2] = minimum_colatitude * constants::degree_to_radians;
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
                    satellite_positions_spherical[p][0] = minimum_radius + ((maximum_radius - minimum_radius) / (n_points_radius - 1)) * h;
                  else
                    satellite_positions_spherical[p][0] = minimum_radius;
                  satellite_positions_spherical[p][2] = std::acos(1. - 2. * s / (n_points_spiral - 1.));
                  satellite_positions_spherical[p][1] = std::fmod((s*golden_angle), 2.*numbers::PI);
                  ++p;
                }
            }
        }
      if (sampling_scheme == list_of_points)
        {
          for (unsigned int p=0; p < n_satellites; ++p)
            {
              if (radius_list.size() == 1)
                satellite_positions_spherical[p][0] = radius_list[0];
              else
                satellite_positions_spherical[p][0] = radius_list[p];
              if (longitude_list[p] < 0)
                satellite_positions_spherical[p][1] = (360 + longitude_list[p]) * constants::degree_to_radians;
              else
                satellite_positions_spherical[p][1] = (longitude_list[p]) * constants::degree_to_radians;
              satellite_positions_spherical[p][2] = (90 - latitude_list[p]) * constants::degree_to_radians;
            }
        }

      // The spherical coordinates are shifted into cartesian to allow simplification
      // in the mathematical equation.
      satellite_positions_cartesian.resize (n_satellites);
      for (unsigned int p=0; p<n_satellites; ++p)
        satellite_positions_cartesian[p]
          = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(satellite_positions_spherical[p]);
    }



    template <int dim>
    std::pair<std::string,std::string>
    GravityPointValues<dim>::execute (TableHandler &)
    {
      AssertThrow(false, ExcNotImplemented());
      return {"", ""};
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
        return {};

      // Up the counter of the number of the file by one, but not in
      // the very first output step. if we run postprocessors on all
      // iterations, only increase file number in the first nonlinear iteration
      const bool increase_file_number = (this->get_nonlinear_iteration() == 0) ||
                                        (!this->get_parameters().run_postprocessors_on_nonlinear_iterations);
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
      const unsigned int degree = this->get_fe().base_element(this->introspection()
                                                              .base_elements.velocities).degree
                                  + quadrature_degree_increase;
      const QGauss<dim> quadrature_formula (degree);
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);


      // Get the value of the universal gravitational constant:
      const double G = aspect::constants::big_g;

      const unsigned int n_quadrature_points_per_cell = quadrature_formula.size();

      const unsigned int n_satellites = satellite_positions_spherical.size();


      // This is the main loop which computes gravity acceleration, potential and
      // gradients at a point located at the spherical coordinate [r, phi, theta].
      // This loop corresponds to the 3 integrals of Newton law:
      std::vector<double>                 local_g_potential (n_satellites);
      std::vector<Tensor<1,dim>>          local_g (n_satellites);
      std::vector<Tensor<1,dim>>          local_g_anomaly (n_satellites);
      std::vector<SymmetricTensor<2,dim>> local_g_gradient (n_satellites);

      std::vector<double> G_times_density_times_JxW (n_quadrature_points_per_cell);
      std::vector<double> G_times_density_anomaly_times_JxW (n_quadrature_points_per_cell);

      MaterialModel::MaterialModelInputs<dim> in(quadrature_formula.size(),
                                                 this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature_formula.size(),
                                                   this->n_compositional_fields());
      in.requested_properties = MaterialModel::MaterialProperties::density;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            // Evaluate the solution at the quadrature points of this cell
            fe_values.reinit (cell);

            in.reinit(fe_values, cell, this->introspection(), this->get_solution());
            this->get_material_model().evaluate(in, out);

            // Pull some computations that are independent of the
            // satellite position out of the following loop. This is
            // because we may have a very large number of satellites,
            // so even just one multiplication that is unnecessarily
            // repeated can be expensive.
            for (unsigned int q = 0; q < n_quadrature_points_per_cell; ++q)
              {
                G_times_density_times_JxW[q]         = G * out.densities[q] *
                                                       fe_values.JxW(q);
                G_times_density_anomaly_times_JxW[q] = G * (out.densities[q]-reference_density) *
                                                       fe_values.JxW(q);
              }


            for (unsigned int p=0; p < n_satellites; ++p)
              {
                const Point<dim> satellite_position = satellite_positions_cartesian[p];

                for (unsigned int q = 0; q < n_quadrature_points_per_cell; ++q)
                  {
                    const Tensor<1,dim> r_vector = satellite_position - fe_values.quadrature_point(q);

                    const double r_squared = r_vector.norm_square();
                    const double r = std::sqrt(r_squared);
                    const double r_cubed = r * r_squared;
                    const double one_over_r_cubed = 1. / r_cubed;
                    const double r_to_the_5 = r_squared * r_cubed;

                    const double G_density_JxW = G_times_density_times_JxW[q];

                    // For gravity acceleration:
                    const double KK = - G_density_JxW * one_over_r_cubed;
                    local_g[p] += KK * r_vector;

                    // For gravity anomalies:
                    const double KK_anomalies = - G_times_density_anomaly_times_JxW[q] * one_over_r_cubed;
                    local_g_anomaly[p] += KK_anomalies * r_vector;

                    // For gravity potential:
                    local_g_potential[p] -= G_density_JxW / r;

                    // For gravity gradient:
                    const double grad_KK = G_density_JxW / r_to_the_5;
                    for (unsigned int e=0; e<dim; ++e)
                      for (unsigned int f=e; f<dim; ++f)
                        local_g_gradient[p][e][f] += grad_KK * (3.0
                                                                * r_vector[e] * r_vector[f]
                                                                - (e==f ? r_squared : 0));
                  }
              }
          }

      // Sum local gravity components over global domain and compute
      // some max and mins. We can directly call Utilities::MPI::sum()
      // for a vector of doubles, but for the other data types we have
      // to be more creative:
      std::vector<double> g_potential(n_satellites);
      Utilities::MPI::sum (make_const_array_view(local_g_potential),
                           this->get_mpi_communicator(),
                           make_array_view(g_potential));

      const auto tensor_sum
        = [] (const auto &v1, const auto &v2)
      {
        AssertDimension (v1.size(), v2.size());

        using element_type = typename std::remove_reference<decltype(v1)>::type::value_type;

        std::vector<element_type> result (v1.size());
        for (unsigned int i=0; i<result.size(); ++i)
          result[i] = v1[i] + v2[i];
        return result;
      };

      const std::vector<Tensor<1,dim>>
      g = Utilities::MPI::all_reduce<decltype(local_g)>
          (local_g,
           this->get_mpi_communicator(),
           tensor_sum);

      const std::vector<Tensor<1,dim>>
      g_anomaly = Utilities::MPI::all_reduce<decltype(local_g_anomaly)>
                  (local_g_anomaly,
                   this->get_mpi_communicator(),
                   tensor_sum);

      const std::vector<SymmetricTensor<2,dim>>
      g_gradient = Utilities::MPI::all_reduce<decltype(local_g_gradient)>
                   (local_g_gradient,
                    this->get_mpi_communicator(),
                    tensor_sum);

      double sum_g = 0;
      double min_g = std::numeric_limits<double>::max();
      double max_g = std::numeric_limits<double>::lowest();
      double sum_g_potential = 0;
      double min_g_potential = std::numeric_limits<double>::max();
      double max_g_potential = std::numeric_limits<double>::lowest();

      for (unsigned int p=0; p < n_satellites; ++p)
        {
          // sum gravity components for all n_satellites:
          sum_g += g[p].norm();
          sum_g_potential += g_potential[p];
          max_g = std::max(g[p].norm(), max_g);
          min_g = std::min(g[p].norm(), min_g);
          max_g_potential = std::max(g_potential[p], max_g_potential);
          min_g_potential = std::min(g_potential[p], min_g_potential);
        }


      // Now also compute theoretical values. These are output only
      // on process zero and, for now, also only computed on process zero.
      std::vector<double> g_theory(n_satellites);
      std::vector<double> g_potential_theory(n_satellites);
      std::vector<SymmetricTensor<2,dim>> g_gradient_theory(n_satellites);

      const unsigned int my_rank
        = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      const unsigned int n_ranks
        = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

      for (unsigned int p=0; p < n_satellites; ++p)
        if (p % n_ranks == my_rank)
          {
            const Point<dim> satellite_position = satellite_positions_cartesian[p];

            // analytical solution to calculate the theoretical gravity and its derivatives
            // from a uniform density model. Can only be used if concentric density profile.
            if (satellite_positions_spherical[p][0] <= model_inner_radius)
              {
                // We are inside the inner radius
                g_theory[p] = 0;
                g_potential_theory[p] = 2.0 * G * numbers::PI * reference_density *
                                        ((model_inner_radius * model_inner_radius) - (model_outer_radius * model_outer_radius));
              }
            else if ((satellite_positions_spherical[p][0] > model_inner_radius)
                     && (satellite_positions_spherical[p][0] < model_outer_radius))
              {
                // We are in the spherical shell
                g_theory[p] = G * numbers::PI * 4./3. * reference_density *
                              (satellite_positions_spherical[p][0] -
                               ((model_inner_radius * model_inner_radius * model_inner_radius)
                                /  (satellite_positions_spherical[p][0] * satellite_positions_spherical[p][0])));
                g_potential_theory[p] = G * numbers::PI * 4./3. * reference_density *
                                        (((satellite_positions_spherical[p][0] * satellite_positions_spherical[p][0])/2.0) +
                                         ((model_inner_radius * model_inner_radius * model_inner_radius) / satellite_positions_spherical[p][0]))
                                        -
                                        G * numbers::PI * 2.0 * reference_density *
                                        (model_outer_radius * model_outer_radius);
              }
            else
              {
                const double common_factor = G * numbers::PI * 4./3. * reference_density
                                             * ((model_outer_radius * model_outer_radius * model_outer_radius) - (model_inner_radius * model_inner_radius * model_inner_radius));
                const double r = satellite_positions_spherical[p][0];

                g_theory[p] = common_factor / (r * r);
                g_potential_theory[p] = - common_factor / r;

                // For the gradient of g, start with the common part of
                // the diagonal elements:
                g_gradient_theory[p][0][0] =
                  g_gradient_theory[p][1][1] =
                    g_gradient_theory[p][2][2] = -1./(r * r * r);

                // Then do the off-diagonal elements:
                for (unsigned int e=0; e<dim; ++e)
                  for (unsigned int f=e; f<dim; ++f)
                    g_gradient_theory[p][e][f] += -(- 3.0 * satellite_position[e] * satellite_position[f])
                                                  /  Utilities::fixed_power<5>(r);
                g_gradient_theory[p] *= common_factor;
              }
          }

      g_theory = Utilities::MPI::all_reduce<decltype(g_theory)>
                 (g_theory,
                  this->get_mpi_communicator(),
                  tensor_sum);

      g_potential_theory = Utilities::MPI::all_reduce<decltype(g_potential_theory)>
                           (g_potential_theory,
                            this->get_mpi_communicator(),
                            tensor_sum);

      g_gradient_theory = Utilities::MPI::all_reduce<decltype(g_gradient_theory)>
                          (g_gradient_theory,
                           this->get_mpi_communicator(),
                           tensor_sum);

      // open the file on rank 0 and write the headers
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // First put all of our output into a string buffer:
          std::ostringstream output;
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

          for (unsigned int p=0; p < n_satellites; ++p)
            {
              const Point<dim> satellite_position = satellite_positions_cartesian[p];

              // write output.
              // g_gradient is here given in eotvos E (1E = 1e-9 per square seconds):
              output << satellite_positions_spherical[p][0] << ' '
                     << satellite_positions_spherical[p][1] * constants::radians_to_degree << ' '
                     << satellite_positions_spherical[p][2] * constants::radians_to_degree << ' '
                     << satellite_position[0] << ' '
                     << satellite_position[1] << ' '
                     << satellite_position[2] << ' '
                     << std::setprecision(precision)
                     << g[p] << ' '
                     << g[p].norm() << ' '
                     << g_theory[p] << ' '
                     << g_potential[p] << ' '
                     << g_potential_theory[p] << ' '
                     << g_anomaly[p] << ' '
                     << g_anomaly[p].norm() << ' '
                     << g_gradient[p][0][0] *1e9 << ' '
                     << g_gradient[p][1][1] *1e9 << ' '
                     << g_gradient[p][2][2] *1e9 << ' '
                     << g_gradient[p][0][1] *1e9 << ' '
                     << g_gradient[p][0][2] *1e9 << ' '
                     << g_gradient[p][1][2] *1e9 << ' '
                     << g_gradient_theory[p][0][0] *1e9 << ' '
                     << g_gradient_theory[p][1][1] *1e9 << ' '
                     << g_gradient_theory[p][2][2] *1e9 << ' '
                     << g_gradient_theory[p][0][1] *1e9 << ' '
                     << g_gradient_theory[p][0][2] *1e9 << ' '
                     << g_gradient_theory[p][1][2] *1e9 << ' '
                     << '\n';
            }

          // Now push the entire buffer into the output file in one swoop fell:
          std::ofstream output_file(filename);
          AssertThrow(output_file,
                      ExcMessage("Unable to open file for writing: " + filename +"."));
          output_file << output.str();
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
                             Patterns::Integer (-1),
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
      const std::string gravity_subdirectory = this->get_output_directory() + "output_gravity/";
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
                       ExcMessage("This postprocessor was only tested for 3d models."));
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
      // Serialize into a stringstream. Put the following into a code
      // block of its own to ensure the destruction of the 'oa'
      // archive triggers a flush() on the stringstream so we can
      // query the completed string below.
      std::ostringstream os;
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

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
