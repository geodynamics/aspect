/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#include <aspect/postprocess/mantle_flux_statistics.h>

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/mpi.h>
#include <deal.II/fe/fe_values.h>

#include <algorithm>
#include <cmath>
#include <queue>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    typename MantleFluxStatistics<dim>::PointDiagnostics
    MantleFluxStatistics<dim>::evaluate_point (const Point<dim> &position,
                                               const Tensor<1,dim> &velocity,
                                               const double temperature,
                                               const double density,
                                               const double thermal_expansivity,
                                               const Tensor<1,dim> &gravity) const
    {
      PointDiagnostics result;
      AssertThrow(gravity.norm() > 0.0,
                  ExcMessage("The mantle flux statistics postprocessor requires "
                             "a nonzero gravity vector."));
      const Tensor<1,dim> outward_direction = -gravity / gravity.norm();

      result.outward_velocity = velocity * outward_direction;
      result.temperature_anomaly =
        temperature - this->get_adiabatic_conditions().temperature(position);

      // Hot outward flow is plume material; cold inward flow is slab material.
      if (result.temperature_anomaly > hot_temperature_threshold &&
          result.outward_velocity > 0.0)
        result.structure = 1;
      else if (result.temperature_anomaly < cold_temperature_threshold &&
               result.outward_velocity < 0.0)
        result.structure = -1;

      if (result.structure != 0)
        {
          // Density and expansivity convert the temperature flux into the
          // anomalous-mass flux commonly called buoyancy flux. Multiplication
          // by gravity then gives the buoyancy force rate.
          result.temperature_flux_density =
            std::abs(result.temperature_anomaly * result.outward_velocity);
          result.buoyancy_mass_flux_density =
            density * thermal_expansivity * result.temperature_flux_density;
          result.buoyancy_force_rate_density =
            gravity.norm() * result.buoyancy_mass_flux_density;
        }

      return result;
    }



    template <int dim>
    double
    MantleFluxStatistics<dim>::surface_distance (const Point<dim> &first,
                                                 const Point<dim> &second) const
    {
      if (this->get_geometry_model().natural_coordinate_system() ==
          Utilities::Coordinates::CoordinateSystem::cartesian)
        {
          const auto first_coordinates =
            this->get_geometry_model().cartesian_to_natural_coordinates(first);
          const auto second_coordinates =
            this->get_geometry_model().cartesian_to_natural_coordinates(second);

          double distance_squared = 0.0;
          for (unsigned int direction = 0; direction < dim-1; ++direction)
            {
              const double difference =
                first_coordinates[direction] - second_coordinates[direction];
              distance_squared += difference * difference;
            }
          return std::sqrt(distance_squared);
        }

      const double first_radius = first.norm();
      const double second_radius = second.norm();
      const double mean_radius = 0.5 * (first_radius + second_radius);
      const double cosine = std::max(-1.0,
                                     std::min(1.0,
                                              (first * second) /
                                              (first_radius * second_radius)));
      return mean_radius * std::acos(cosine);
    }



    template <int dim>
    std::vector<std::vector<Point<dim>>>
    MantleFluxStatistics<dim>::find_clusters (const std::vector<Point<dim>> &points) const
    {
      std::vector<std::vector<Point<dim>>> clusters;
      std::vector<bool> visited(points.size(), false);

      for (unsigned int first_point = 0; first_point < points.size(); ++first_point)
        {
          if (visited[first_point])
            continue;

          std::queue<unsigned int> points_to_visit;
          std::vector<Point<dim>> cluster;
          points_to_visit.push(first_point);

          while (!points_to_visit.empty())
            {
              const unsigned int current_point = points_to_visit.front();
              points_to_visit.pop();

              if (visited[current_point])
                continue;

              visited[current_point] = true;
              cluster.push_back(points[current_point]);

              for (unsigned int neighbor = 0; neighbor < points.size(); ++neighbor)
                if (!visited[neighbor] &&
                    surface_distance(points[current_point], points[neighbor]) <= maximum_point_spacing)
                  points_to_visit.push(neighbor);
            }

          if (cluster.size() >= minimum_points_per_structure)
            clusters.push_back(cluster);
        }

      return clusters;
    }



    template <int dim>
    double
    MantleFluxStatistics<dim>::cluster_length (const std::vector<Point<dim>> &points) const
    {
      double length = 0.0;
      for (unsigned int first = 0; first < points.size(); ++first)
        for (unsigned int second = first + 1; second < points.size(); ++second)
          length = std::max(length, surface_distance(points[first], points[second]));
      return length;
    }



    template <int dim>
    double
    MantleFluxStatistics<dim>::cluster_radius (const std::vector<Point<dim>> &points) const
    {
      Point<dim> center;
      for (const Point<dim> &point : points)
        center += point;
      center /= points.size();

      if (this->get_geometry_model().natural_coordinate_system() ==
          Utilities::Coordinates::CoordinateSystem::spherical)
        {
          Point<dim> direction_center;
          double mean_shell_radius = 0.0;
          for (const Point<dim> &point : points)
            {
              direction_center += point / point.norm();
              mean_shell_radius += point.norm();
            }

          mean_shell_radius /= points.size();
          if (direction_center.norm() == 0.0)
            return 0.0;
          center = direction_center * (mean_shell_radius / direction_center.norm());
        }

      double radius = 0.0;
      for (const Point<dim> &point : points)
        radius += surface_distance(center, point);
      return radius / points.size();
    }



    template <int dim>
    void
    MantleFluxStatistics<dim>::write_cluster_file (const std::string &content)
    {
      if (std::isnan(last_output_time))
        last_output_time = this->get_time() - cluster_file_interval;

      if (!write_cluster_files ||
          ((this->get_time() < last_output_time + cluster_file_interval) &&
           (this->get_timestep_number() != 0)))
        return;

      std::string filename = this->get_output_directory() +
                             "mantle_flux_clusters." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);

      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename += "." + Utilities::int_to_string(this->get_nonlinear_iteration(), 4);

      Utilities::collect_and_write_file_content(filename,
                                                content,
                                                this->get_mpi_communicator());

      if (cluster_file_interval > 0.0)
        {
          const double adjustment = 1.0 + 2.0 * std::numeric_limits<double>::epsilon();
          last_output_time +=
            std::floor((this->get_time() - last_output_time) /
                       cluster_file_interval * adjustment) *
            cluster_file_interval / adjustment;
        }
    }



    template <int dim>
    std::pair<std::string,std::string>
    MantleFluxStatistics<dim>::execute (TableHandler &statistics)
    {
      const auto coordinate_system =
        this->get_geometry_model().natural_coordinate_system();
      AssertThrow(coordinate_system == Utilities::Coordinates::CoordinateSystem::spherical ||
                  coordinate_system == Utilities::Coordinates::CoordinateSystem::cartesian,
                  ExcMessage("The mantle flux statistics postprocessor supports "
                             "spherical and Cartesian geometry models."));

      const MPI_Comm communicator = this->get_mpi_communicator();
      const unsigned int rank = Utilities::MPI::this_mpi_process(communicator);
      const Quadrature<dim> &quadrature = this->introspection().quadratures.velocities;
      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature,
                              update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> material_inputs(quadrature.size(),
                                                               this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> material_outputs(quadrature.size(),
                                                                 this->n_compositional_fields());
      material_inputs.requested_properties =
        MaterialModel::MaterialProperties::density |
        MaterialModel::MaterialProperties::thermal_expansion_coefficient;

      std::ostringstream screen_output;
      std::ostringstream cluster_output;

      if (rank == 0)
        {
          cluster_output << "# ";
          for (unsigned int direction = 0; direction < dim; ++direction)
            cluster_output << static_cast<char>('x' + direction) << ' ';
          cluster_output << "depth_km type cluster_id\n";
        }

      unsigned int next_cluster_id = 0;

      enum FluxComponent
        {
          plume_volume_flux,
          slab_volume_flux,
          plume_temperature_flux,
          slab_temperature_flux,
          plume_buoyancy_mass_flux,
          slab_buoyancy_mass_flux,
          plume_buoyancy_force_rate,
          slab_buoyancy_force_rate,
          n_flux_components
        };

      for (const double requested_depth : measurement_depths)
        {
          std::vector<double> local_fluxes(n_flux_components, 0.0);
          std::vector<Point<dim>> local_plume_points;
          std::vector<Point<dim>> local_slab_points;

          // Sample a thin layer around the requested depth. Dividing its
          // volume by the layer thickness approximates a surface integral.
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);
                material_inputs.reinit(fe_values, cell, this->introspection(), this->get_solution());
                this->get_material_model().evaluate(material_inputs, material_outputs);

                bool cell_contains_plume = false;
                bool cell_contains_slab = false;

                for (unsigned int q = 0; q < quadrature.size(); ++q)
                  {
                    const Point<dim> &position = fe_values.quadrature_point(q);
                    if (std::abs(this->get_geometry_model().depth(position) - requested_depth) >
                        measurement_half_thickness)
                      continue;

                    const Tensor<1,dim> gravity =
                      this->get_gravity_model().gravity_vector(position);
                    const PointDiagnostics point =
                      evaluate_point(position,
                                     material_inputs.velocity[q],
                                     material_inputs.temperature[q],
                                     material_outputs.densities[q],
                                     material_outputs.thermal_expansion_coefficients[q],
                                     gravity);
                    const double surface_measure =
                      fe_values.JxW(q) / (2.0 * measurement_half_thickness);

                    if (point.structure == 1)
                      {
                        local_fluxes[plume_volume_flux] +=
                          point.outward_velocity * surface_measure;
                        local_fluxes[plume_temperature_flux] +=
                          point.temperature_flux_density * surface_measure;
                        local_fluxes[plume_buoyancy_mass_flux] +=
                          point.buoyancy_mass_flux_density * surface_measure;
                        local_fluxes[plume_buoyancy_force_rate] +=
                          point.buoyancy_force_rate_density * surface_measure;
                        cell_contains_plume = true;
                      }
                    else if (point.structure == -1)
                      {
                        local_fluxes[slab_volume_flux] +=
                          -point.outward_velocity * surface_measure;
                        local_fluxes[slab_temperature_flux] +=
                          point.temperature_flux_density * surface_measure;
                        local_fluxes[slab_buoyancy_mass_flux] +=
                          point.buoyancy_mass_flux_density * surface_measure;
                        local_fluxes[slab_buoyancy_force_rate] +=
                          point.buoyancy_force_rate_density * surface_measure;
                        cell_contains_slab = true;
                      }
                  }

                if (cell_contains_plume)
                  local_plume_points.push_back(cell->center());
                if (cell_contains_slab)
                  local_slab_points.push_back(cell->center());
              }

          // Collect all integrated fluxes in one communication step.
          std::vector<double> global_fluxes(n_flux_components);
          Utilities::MPI::sum(local_fluxes, communicator, global_fluxes);

          // Gather detected cell centers so that connected plumes and slabs
          // can be counted consistently across MPI boundaries.
          const std::vector<std::vector<Point<dim>>> gathered_plume_points =
            Utilities::MPI::gather(communicator, local_plume_points);
          const std::vector<std::vector<Point<dim>>> gathered_slab_points =
            Utilities::MPI::gather(communicator, local_slab_points);

          unsigned int plume_count = 0;
          unsigned int slab_count = 0;
          double total_slab_length = 0.0;
          double mean_slab_length = 0.0;
          double mean_plume_radius = 0.0;

          if (rank == 0)
            {
              std::vector<Point<dim>> plume_points;
              std::vector<Point<dim>> slab_points;
              for (const auto &points : gathered_plume_points)
                plume_points.insert(plume_points.end(), points.begin(), points.end());
              for (const auto &points : gathered_slab_points)
                slab_points.insert(slab_points.end(), points.begin(), points.end());

              // Nearby points form one structure. Slabs shorter than the
              // requested minimum length are left out of the reported count.
              const auto plume_clusters = find_clusters(plume_points);
              const auto possible_slab_clusters = find_clusters(slab_points);
              std::vector<std::vector<Point<dim>>> slab_clusters;
              for (const auto &cluster : possible_slab_clusters)
                if (cluster_length(cluster) >= minimum_slab_length)
                  slab_clusters.push_back(cluster);

              plume_count = plume_clusters.size();
              slab_count = slab_clusters.size();

              for (const auto &cluster : plume_clusters)
                {
                  mean_plume_radius += cluster_radius(cluster);
                  const unsigned int cluster_id = next_cluster_id++;
                  for (const Point<dim> &point : cluster)
                    {
                      for (unsigned int direction = 0; direction < dim; ++direction)
                        cluster_output << point[direction] << ' ';
                      cluster_output << requested_depth / 1000.0 << " plume " << cluster_id << '\n';
                    }
                }

              for (const auto &cluster : slab_clusters)
                {
                  total_slab_length += cluster_length(cluster);
                  const unsigned int cluster_id = next_cluster_id++;
                  for (const Point<dim> &point : cluster)
                    {
                      for (unsigned int direction = 0; direction < dim; ++direction)
                        cluster_output << point[direction] << ' ';
                      cluster_output << requested_depth / 1000.0 << " slab " << cluster_id << '\n';
                    }
                }

              if (plume_count > 0)
                mean_plume_radius /= plume_count;
              if (slab_count > 0)
                mean_slab_length = total_slab_length / slab_count;
            }

          plume_count = Utilities::MPI::broadcast(communicator, plume_count, 0);
          slab_count = Utilities::MPI::broadcast(communicator, slab_count, 0);
          total_slab_length = Utilities::MPI::broadcast(communicator, total_slab_length, 0);
          mean_slab_length = Utilities::MPI::broadcast(communicator, mean_slab_length, 0);
          mean_plume_radius = Utilities::MPI::broadcast(communicator, mean_plume_radius, 0);

          // A 2D flux is per unit length in the missing direction; a 3D flux
          // is a volume per time.
          const bool output_in_years = this->convert_output_to_years();
          const double time_conversion = output_in_years ? year_in_seconds : 1.0;
          const std::string time_unit = output_in_years ? "yr" : "s";
          const double volume_conversion = time_conversion / std::pow(1000.0, dim);
          const std::string volume_unit =
            (dim == 3 ? "km^3/" : "km^2/") + time_unit;
          const std::string temperature_flux_unit =
            (dim == 3 ? "K km^3/" : "K km^2/") + time_unit;
          const std::string buoyancy_mass_flux_unit =
            (dim == 3 ? "kg/" + time_unit : "kg/(m " + time_unit + ")");
          const std::string buoyancy_unit =
            (dim == 3 ? "N/" + time_unit : "N/(m " + time_unit + ")");
          const std::string depth_label = Utilities::to_string(requested_depth / 1000.0) + " km";

          const auto add_scientific_value = [&statistics](const std::string &name, const double value)
          {
            statistics.add_value(name, value);
            statistics.set_precision(name, 8);
            statistics.set_scientific(name, true);
          };

          add_scientific_value("Plume volume flux at " + depth_label + " (" + volume_unit + ")",
                               global_fluxes[plume_volume_flux] * volume_conversion);
          add_scientific_value("Slab volume flux at " + depth_label + " (" + volume_unit + ")",
                               global_fluxes[slab_volume_flux] * volume_conversion);
          add_scientific_value("Plume temperature anomaly flux at " + depth_label + " (" + temperature_flux_unit + ")",
                               global_fluxes[plume_temperature_flux] * volume_conversion);
          add_scientific_value("Slab temperature anomaly flux at " + depth_label + " (" + temperature_flux_unit + ")",
                               global_fluxes[slab_temperature_flux] * volume_conversion);
          add_scientific_value("Plume thermal buoyancy mass flux at " + depth_label + " (" + buoyancy_mass_flux_unit + ")",
                               global_fluxes[plume_buoyancy_mass_flux] * time_conversion);
          add_scientific_value("Slab thermal buoyancy mass flux at " + depth_label + " (" + buoyancy_mass_flux_unit + ")",
                               global_fluxes[slab_buoyancy_mass_flux] * time_conversion);
          add_scientific_value("Plume thermal buoyancy force rate at " + depth_label + " (" + buoyancy_unit + ")",
                               global_fluxes[plume_buoyancy_force_rate] * time_conversion);
          add_scientific_value("Slab thermal buoyancy force rate at " + depth_label + " (" + buoyancy_unit + ")",
                               global_fluxes[slab_buoyancy_force_rate] * time_conversion);
          statistics.add_value("Number of plumes at " + depth_label, plume_count);
          statistics.add_value("Number of slabs at " + depth_label, slab_count);
          add_scientific_value("Total slab length at " + depth_label + " (km)", total_slab_length / 1000.0);
          add_scientific_value("Mean slab length at " + depth_label + " (km)", mean_slab_length / 1000.0);
          add_scientific_value("Mean plume radius at " + depth_label + " (km)", mean_plume_radius / 1000.0);

          screen_output << "Depth " << depth_label
                        << ": plume flux="
                        << global_fluxes[plume_volume_flux] * volume_conversion << ' ' << volume_unit
                        << ", slab flux="
                        << global_fluxes[slab_volume_flux] * volume_conversion << ' ' << volume_unit
                        << ", plume temperature flux="
                        << global_fluxes[plume_temperature_flux] * volume_conversion << ' '
                        << temperature_flux_unit
                        << ", slab temperature flux="
                        << global_fluxes[slab_temperature_flux] * volume_conversion << ' '
                        << temperature_flux_unit
                        << ", plume buoyancy mass flux="
                        << global_fluxes[plume_buoyancy_mass_flux] * time_conversion << ' '
                        << buoyancy_mass_flux_unit
                        << ", slab buoyancy mass flux="
                        << global_fluxes[slab_buoyancy_mass_flux] * time_conversion << ' '
                        << buoyancy_mass_flux_unit
                        << ", plume buoyancy force rate="
                        << global_fluxes[plume_buoyancy_force_rate] * time_conversion << ' ' << buoyancy_unit
                        << ", slab buoyancy force rate="
                        << global_fluxes[slab_buoyancy_force_rate] * time_conversion << ' ' << buoyancy_unit
                        << ", plumes=" << plume_count
                        << ", slabs=" << slab_count << ". ";
        }

      write_cluster_file(cluster_output.str());
      return {"Mantle flux statistics:", screen_output.str()};
    }



    template <int dim>
    void
    MantleFluxStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Mantle flux statistics");
        {
          prm.declare_entry("Measurement depths", "440e3",
                            Patterns::List(Patterns::Double(0.0)),
                            "Depths below the surface where plume and slab fluxes are measured. Units: m.");
          prm.declare_entry("Measurement layer half thickness", "20000",
                            Patterns::Double(0.0),
                            "Half the thickness of the layer sampled around each measurement depth. Units: m.");
          prm.declare_entry("Cold temperature anomaly threshold", "-200",
                            Patterns::Double(),
                            "A point colder than this value and moving inward is treated as slab material. Units: K.");
          prm.declare_entry("Hot temperature anomaly threshold", "200",
                            Patterns::Double(),
                            "A point hotter than this value and moving outward is "
                            "treated as plume material. Units: K.");
          prm.declare_entry("Maximum point spacing", "80000",
                            Patterns::Double(0.0),
                            "Largest surface distance between neighboring cell centers "
                            "that belong to the same plume or slab. Units: m.");
          prm.declare_entry("Minimum points per structure", "5",
                            Patterns::Integer(1),
                            "Smallest number of detected cell centers needed to count a plume or slab.");
          prm.declare_entry("Minimum slab length", "0",
                            Patterns::Double(0.0),
                            "Do not count cold structures shorter than this surface distance. Units: m.");
          prm.declare_entry("Write cluster files", "false",
                            Patterns::Bool(),
                            "Write the cell centers assigned to each plume and slab "
                            "to mantle_flux_clusters.NNNNN files.");
          prm.declare_entry("Time between cluster files", "0",
                            Patterns::Double(0.0),
                            "Time between cluster files. Zero writes a file every time "
                            "the postprocessor runs. Units: years when output uses years; "
                            "seconds otherwise.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MantleFluxStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Mantle flux statistics");
        {
          measurement_depths = Utilities::string_to_double(
            Utilities::split_string_list(prm.get("Measurement depths")));
          measurement_half_thickness = prm.get_double("Measurement layer half thickness");
          cold_temperature_threshold = prm.get_double("Cold temperature anomaly threshold");
          hot_temperature_threshold = prm.get_double("Hot temperature anomaly threshold");
          maximum_point_spacing = prm.get_double("Maximum point spacing");
          minimum_points_per_structure = prm.get_integer("Minimum points per structure");
          minimum_slab_length = prm.get_double("Minimum slab length");
          write_cluster_files = prm.get_bool("Write cluster files");
          cluster_file_interval = prm.get_double("Time between cluster files");
          if (this->convert_output_to_years())
            cluster_file_interval *= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      AssertThrow(measurement_half_thickness > 0.0,
                  ExcMessage("Measurement layer half thickness must be greater than zero."));
      AssertThrow(cold_temperature_threshold < 0.0,
                  ExcMessage("Cold temperature anomaly threshold must be negative."));
      AssertThrow(hot_temperature_threshold > 0.0,
                  ExcMessage("Hot temperature anomaly threshold must be positive."));
      AssertThrow(maximum_point_spacing > 0.0,
                  ExcMessage("Maximum point spacing must be greater than zero."));
    }
  }
}



namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MantleFluxStatistics,
                                  "mantle flux statistics",
                                  "Measures hot outward flow (plumes) and cold inward flow (slabs) "
                                  "across layers at selected depths in spherical and Cartesian geometries. It reports "
                                  "volume flux, temperature-anomaly flux, thermal-buoyancy mass flux, "
                                  "thermal-buoyancy force rate, the "
                                  "number of detected structures, slab length, and plume radius. In "
                                  "two-dimensional models, fluxes are reported per unit length in the "
                                  "missing third direction. Rate units follow the global choice to "
                                  "use years or seconds in output.")
  }
}
