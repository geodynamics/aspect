/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE. If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/continental_fragmentation_statistics.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>

namespace aspect
{
  namespace Postprocess
  {
    // 1. FaceRecord and BlockSummary store the measured properties of each
    // surface face and connected continental block.
    namespace
    {
      // 2. Provide geometric operations for matching vertices, measuring edge
      // lengths, and reporting surface positions as longitude and latitude.
      template <int dim>
      double point_distance(const Point<dim> &first_point,
                            const Point<dim> &second_point)
      {
        return first_point.distance(second_point);
      }


      template <int dim>
      std::string point_key(const Point<dim> &point,
                            const double tolerance = 1e-10)
      {
        // Match repeated copies of the same surface vertex when reconstructing
        // the global continent geometry.
        std::ostringstream out;
        out << std::fixed << std::setprecision(0);

        for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
          {
            const double quantized_coordinate =
              std::round(point[coordinate] / tolerance);
            out << quantized_coordinate;
            if (coordinate + 1 < dim)
              out << "_";
          }

        return out.str();
      }


      template <int dim>
      std::pair<double,double> cartesian_to_lon_lat_deg(const Point<dim> &point)
      {
        static_assert(dim == 3, "lon/lat conversion only implemented for dim=3.");

        const double radius = point.norm();
        if (radius <= 0.0)
          return std::make_pair(0.0, 0.0);

        const double longitude =
          std::atan2(point[1], point[0]) * 180.0 / numbers::PI;
        const double latitude =
          std::asin(point[2] / radius) * 180.0 / numbers::PI;
        return std::make_pair(longitude, latitude);
      }

      // 3. Use a union-find structure to group continental surface faces that
      // share a vertex into connected blocks at the current timestep.
      class DisjointSet
      {
        public:
          explicit DisjointSet(const unsigned int size)
            : parent_indices(size), tree_depths(size, 0)
          {
            for (unsigned int index = 0; index < size; ++index)
              parent_indices[index] = index;
          }

          unsigned int find_root(const unsigned int index)
          {
            if (parent_indices[index] != index)
              parent_indices[index] = find_root(parent_indices[index]);
            return parent_indices[index];
          }

          void merge(const unsigned int first_index,
                     const unsigned int second_index)
          {
            unsigned int first_root = find_root(first_index);
            unsigned int second_root = find_root(second_index);

            if (first_root == second_root)
              return;

            if (tree_depths[first_root] < tree_depths[second_root])
              std::swap(first_root, second_root);

            parent_indices[second_root] = first_root;

            if (tree_depths[first_root] == tree_depths[second_root])
              ++tree_depths[first_root];
          }

        private:
          std::vector<unsigned int> parent_indices;
          std::vector<unsigned int> tree_depths;
      };


      using EdgeKey = std::pair<std::string,std::string>;

      struct EdgeRecord
      {
        std::vector<unsigned int> incident_face_indices;
        double length = 0.0;
      };


    }



    template <int dim>
    void
    ContinentalFragmentationStatistics<dim>::declare_parameters(ParameterHandler &prm)
    {
      // 4. Define the continental composition fields, classification threshold,
      // minimum block area, and optional diagnostic outputs.
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Continental fragmentation statistics");
        {
          prm.declare_entry("Continent field names",
                            "",
                            Patterns::List(Patterns::Anything()),
                            "Comma-separated list of compositional field names "
                            "that should be interpreted as continental material.");

          prm.declare_entry("Continent threshold",
                            "0.5",
                            Patterns::Double(0.0),
                            "A top-surface face is considered continental if the "
                            "sum of the selected compositional fields exceeds this threshold.");

          prm.declare_entry("Minimum block area",
                            "0.0",
                            Patterns::Double(0.0),
                            "Ignore connected continental blocks with area smaller "
                            "than this value [m^2].");

          prm.declare_entry("Output verbose screen line",
                            "true",
                            Patterns::Bool(),
                            "Whether to print a compact summary line to the screen.");

          prm.declare_entry("Write surface map",
                            "true",
                            Patterns::Bool(),
                            "Write a per-face diagnostic file for the top surface.");

          prm.declare_entry("Write block summary",
                            "true",
                            Patterns::Bool(),
                            "Write a per-block diagnostic summary file.");

          prm.declare_entry("Output file prefix",
                            "continental_fragmentation_statistics",
                            Patterns::Anything(),
                            "Prefix used for the diagnostic output file names.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ContinentalFragmentationStatistics<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Continental fragmentation statistics");
        {
          continent_field_names =
            Utilities::split_string_list(prm.get("Continent field names"));

          continent_threshold = prm.get_double("Continent threshold");
          minimum_block_area = prm.get_double("Minimum block area");
          output_verbose_screen_line = prm.get_bool("Output verbose screen line");

          write_surface_map = prm.get_bool("Write surface map");
          write_block_summary = prm.get_bool("Write block summary");
          output_file_prefix = prm.get("Output file prefix");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      resolve_field_indices();
    }



    template <int dim>
    void
    ContinentalFragmentationStatistics<dim>::resolve_field_indices()
    {
      continent_field_indices.clear();

      const std::vector<std::string> &all_names =
        this->introspection().chemical_composition_field_names();

      for (const std::string &requested_name : continent_field_names)
        {
          bool found = false;

          for (unsigned int field_index = 0;
               field_index < all_names.size();
               ++field_index)
            {
              if (all_names[field_index] == requested_name)
                {
                  continent_field_indices.push_back(field_index);
                  found = true;
                  break;
                }
            }

          AssertThrow(found,
                      ExcMessage("Continental fragmentation statistics: field <" +
                                 requested_name +
                                 "> was not found among the compositional fields."));
        }

      AssertThrow(!continent_field_indices.empty(),
                  ExcMessage("Continental fragmentation statistics: at least one continent field "
                             "must be provided."));
    }



    template <int dim>
    bool
    ContinentalFragmentationStatistics<dim>::is_top_boundary_face(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      const unsigned int face_no) const
    {
      // 5. Restrict the analysis to the physical top surface of the model.
      if (!cell->face(face_no)->at_boundary())
        return false;

      const types::boundary_id top_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      return (cell->face(face_no)->boundary_id() == top_id);
    }


    template <int dim>
    void
    ContinentalFragmentationStatistics<dim>::write_surface_file(const std::vector<FaceRecord> &all_faces) const
    {
      std::ostringstream filename;
      filename << this->get_output_directory()
               << output_file_prefix
               << "_surface."
               << std::setw(5) << std::setfill('0')
               << this->get_timestep_number();

      std::ofstream out(filename.str().c_str());
      AssertThrow(out,
                  ExcMessage("Could not open continent surface output file: " + filename.str()));

      out << "# 1:x 2:y ";
      if constexpr (dim == 3)
        out << "3:z 4:radius 5:lon_deg 6:lat_deg 7:area 8:continent_value 9:is_continent 10:speed 11:block_id\n";
      else
        out << "3:area 4:continent_value 5:is_continent 6:speed 7:block_id\n";

      out << std::setprecision(16);

      for (const FaceRecord &surface_face : all_faces)
        {
          if constexpr (dim == 3)
            {
              const auto longitude_latitude =
                cartesian_to_lon_lat_deg(surface_face.center);

              out << surface_face.center[0] << " "
                  << surface_face.center[1] << " "
                  << surface_face.center[2] << " "
                  << surface_face.center.norm() << " "
                  << longitude_latitude.first << " "
                  << longitude_latitude.second << " "
                  << surface_face.area << " "
                  << surface_face.continent_value << " "
                  << (surface_face.is_continent ? 1 : 0) << " "
                  << surface_face.speed << " "
                  << surface_face.block_id << "\n";
            }
          else
            {
              out << surface_face.center[0] << " "
                  << surface_face.center[1] << " "
                  << surface_face.area << " "
                  << surface_face.continent_value << " "
                  << (surface_face.is_continent ? 1 : 0) << " "
                  << surface_face.speed << " "
                  << surface_face.block_id << "\n";
            }
        }
    }



    template <int dim>
    void
    ContinentalFragmentationStatistics<dim>::write_block_file(
      const std::vector<BlockSummary> &blocks,
      const std::vector<Point<dim>> &block_centroids) const
    {
      std::ostringstream filename;
      filename << this->get_output_directory()
               << output_file_prefix
               << "_blocks."
               << std::setw(5) << std::setfill('0')
               << this->get_timestep_number();

      std::ofstream out(filename.str().c_str());
      AssertThrow(out,
                  ExcMessage("Could not open continent block output file: " + filename.str()));

      out << "# 1:block_id ";
      if constexpr (dim == 3)
        out << "2:centroid_x 3:centroid_y 4:centroid_z 5:radius 6:lon_deg 7:lat_deg 8:area 9:mean_speed 10:n_faces\n";
      else
        out << "2:centroid_x 3:centroid_y 4:area 5:mean_speed 6:n_faces\n";

      out << std::setprecision(16);

      for (unsigned int block_index = 0;
           block_index < blocks.size();
           ++block_index)
        {
          const double mean_speed =
            (blocks[block_index].area > 0.0
             ? blocks[block_index].speed_area_integral / blocks[block_index].area
             : 0.0);

          if constexpr (dim == 3)
            {
              const auto longitude_latitude =
                cartesian_to_lon_lat_deg(block_centroids[block_index]);

              out << blocks[block_index].id << " "
                  << block_centroids[block_index][0] << " "
                  << block_centroids[block_index][1] << " "
                  << block_centroids[block_index][2] << " "
                  << block_centroids[block_index].norm() << " "
                  << longitude_latitude.first << " "
                  << longitude_latitude.second << " "
                  << blocks[block_index].area << " "
                  << mean_speed << " "
                  << blocks[block_index].n_faces << "\n";
            }
          else
            {
              out << blocks[block_index].id << " "
                  << block_centroids[block_index][0] << " "
                  << block_centroids[block_index][1] << " "
                  << blocks[block_index].area << " "
                  << mean_speed << " "
                  << blocks[block_index].n_faces << "\n";
            }
        }
    }



    template <int dim>
    std::pair<std::string,std::string>
    ContinentalFragmentationStatistics<dim>::execute(TableHandler &statistics)
    {
      resolve_field_indices();

      const QGauss<dim-1> face_quadrature(this->get_fe().base_element(0).degree + 1);

      FEFaceValues<dim> face_values(this->get_mapping(),
                                    this->get_fe(),
                                    face_quadrature,
                                    update_values |
                                    update_quadrature_points |
                                    update_JxW_values);

      std::vector<Tensor<1,dim>> velocity_values(face_quadrature.size());

      std::vector<std::vector<double>> composition_values_per_field(
        continent_field_indices.size(),
        std::vector<double>(face_quadrature.size()));

      std::vector<FaceRecord> local_faces;

      const auto &velocity_extractor = this->introspection().extractors.velocities;

      // 6. Sample the top surface to measure each face's area, mean continent
      // fraction, and mean drift speed.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (!cell->is_locally_owned())
            continue;

          for (const unsigned int face_number : cell->face_indices())
            {
              if (!is_top_boundary_face(cell, face_number))
                continue;

              face_values.reinit(cell, face_number);

              face_values[velocity_extractor].get_function_values(this->get_solution(),
                                                                  velocity_values);

              for (unsigned int selected_field = 0;
                   selected_field < continent_field_indices.size();
                   ++selected_field)
                {
                  const FEValuesExtractors::Scalar composition_extractor =
                    this->introspection().extractors.compositional_fields[
                      continent_field_indices[selected_field]];

                  face_values[composition_extractor].get_function_values(
                    this->get_solution(),
                    composition_values_per_field[selected_field]);
                }

              double face_area = 0.0;
              double face_speed_integral = 0.0;
              double face_continent_integral = 0.0;

              for (unsigned int quadrature_point = 0;
                   quadrature_point < face_quadrature.size();
                   ++quadrature_point)
                {
                  double summed_continent_composition = 0.0;
                  for (unsigned int selected_field = 0;
                       selected_field < continent_field_indices.size();
                       ++selected_field)
                    summed_continent_composition +=
                      composition_values_per_field[selected_field][quadrature_point];

                  const double quadrature_weight = face_values.JxW(quadrature_point);
                  face_area += quadrature_weight;
                  face_speed_integral +=
                    velocity_values[quadrature_point].norm() * quadrature_weight;
                  face_continent_integral +=
                    summed_continent_composition * quadrature_weight;
                }

              FaceRecord face_record;
              face_record.center = cell->face(face_number)->center();
              face_record.area = face_area;
              face_record.continent_value =
                face_continent_integral / std::max(face_area, 1e-30);
              face_record.speed =
                face_speed_integral / std::max(face_area, 1e-30 * year_in_seconds);

              face_record.is_continent =
                (face_record.continent_value > continent_threshold);

              face_record.vertices.reserve(GeometryInfo<dim>::vertices_per_face);
              for (unsigned int vertex_index = 0;
                   vertex_index < GeometryInfo<dim>::vertices_per_face;
                   ++vertex_index)
                face_record.vertices.push_back(
                  cell->face(face_number)->vertex(vertex_index));

              local_faces.push_back(face_record);
            }
        }

      std::vector<double> packed_local;
      packed_local.reserve(local_faces.size() *
                           (4 + dim + dim * GeometryInfo<dim>::vertices_per_face));

      for (const FaceRecord &face_record : local_faces)
        {
          packed_local.push_back(face_record.area);
          packed_local.push_back(face_record.continent_value);
          packed_local.push_back(face_record.speed);
          packed_local.push_back(face_record.is_continent ? 1.0 : 0.0);

          for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
            packed_local.push_back(face_record.center[coordinate]);

          for (const Point<dim> &vertex : face_record.vertices)
            for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
              packed_local.push_back(vertex[coordinate]);
        }

      // 7. Assemble one global surface map so continent blocks can cross MPI
      // partition boundaries.
      const std::vector<std::vector<double>> gathered =
        Utilities::MPI::gather(this->get_mpi_communicator(), packed_local, 0);

      double global_cont_area = 0.0;
      double global_cont_speed_integral = 0.0;
      double global_perimeter = 0.0;
      double global_largest_block = 0.0;
      unsigned int global_n_blocks = 0;
      double fragmentation_index = 0.0;
      double perimeter_to_area_ratio = 0.0;
      double normalized_perimeter_fragmentation = 0.0;

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Reconstruct the complete top-surface map from all MPI processes.
          std::vector<FaceRecord> all_faces;

          // Each packed face contains its area, composition, velocity,
          // classification, center, and vertex coordinates.

          for (const auto &buffer : gathered)
            {
              const unsigned int stride =
                4 + dim + dim * GeometryInfo<dim>::vertices_per_face;

              AssertThrow(buffer.size() % stride == 0,
                          ExcMessage("Continental fragmentation statistics: invalid gathered buffer size."));

              unsigned int buffer_offset = 0;
              const unsigned int number_of_buffered_faces = buffer.size() / stride;

              for (unsigned int face_index = 0;
                   face_index < number_of_buffered_faces;
                   ++face_index)
                {
                  FaceRecord face_record;
                  face_record.area = buffer[buffer_offset++];
                  face_record.continent_value = buffer[buffer_offset++];
                  face_record.speed = buffer[buffer_offset++];
                  face_record.is_continent = (buffer[buffer_offset++] > 0.5);

                  for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
                    face_record.center[coordinate] = buffer[buffer_offset++];

                  face_record.vertices.resize(GeometryInfo<dim>::vertices_per_face);
                  for (unsigned int vertex_index = 0;
                       vertex_index < GeometryInfo<dim>::vertices_per_face;
                       ++vertex_index)
                    for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
                      face_record.vertices[vertex_index][coordinate] =
                        buffer[buffer_offset++];

                  all_faces.push_back(face_record);
                }
            }

          const unsigned int number_of_faces = all_faces.size();

          // 8. Integrate total continental area and area-weighted drift speed.
          for (const FaceRecord &surface_face : all_faces)
            {
              if (surface_face.is_continent)
                {
                  global_cont_area += surface_face.area;
                  global_cont_speed_integral += surface_face.speed * surface_face.area;
                }
            }

          std::map<std::string, std::vector<unsigned int>> vertex_to_faces;
          for (unsigned int face_index = 0; face_index < number_of_faces; ++face_index)
            for (const Point<dim> &vertex : all_faces[face_index].vertices)
              vertex_to_faces[point_key(vertex)].push_back(face_index);

          // 9. Identify continent blocks as connected sets of surface faces
          // whose mean continent fraction exceeds the threshold.
          // Vertex A → faces 10, 11, 15
          // Vertex B → faces 11, 12
          DisjointSet connected_faces(number_of_faces);

          for (const auto &vertex_entry : vertex_to_faces)
            {
              const std::vector<unsigned int> &incident_face_indices = vertex_entry.second;

              for (unsigned int first_position = 0;
                   first_position < incident_face_indices.size();
                   ++first_position)
                for (unsigned int second_position = first_position + 1;
                     second_position < incident_face_indices.size();
                     ++second_position)
                  {
                    const unsigned int first_face_index =
                      incident_face_indices[first_position];
                    const unsigned int second_face_index =
                      incident_face_indices[second_position];

                    if (all_faces[first_face_index].is_continent &&
                        all_faces[second_face_index].is_continent)
                      connected_faces.merge(first_face_index, second_face_index);
                  }
            }

          // 10. Measure each connected block and assign IDs only to blocks above
          // the user-defined minimum block area. Smaller blocks remain part of
          // the total continental area and perimeter.
          std::map<unsigned int, double> component_areas;
          for (unsigned int face_index = 0; face_index < number_of_faces; ++face_index)
            if (all_faces[face_index].is_continent)
              {
                const unsigned int component_root =
                  connected_faces.find_root(face_index);
                component_areas[component_root] += all_faces[face_index].area;
              }

          std::map<unsigned int, unsigned int> component_root_to_block_id;
          unsigned int next_block_id = 0;

          for (const auto &component : component_areas)
            {
              const unsigned int component_root = component.first;
              const double component_area = component.second;

              if (component_area >= minimum_block_area)
                {
                  component_root_to_block_id[component_root] = next_block_id++;
                  ++global_n_blocks;
                  global_largest_block = std::max(global_largest_block, component_area);
                }
            }

          for (unsigned int face_index = 0; face_index < number_of_faces; ++face_index)
            {
              if (!all_faces[face_index].is_continent)
                {
                  all_faces[face_index].block_id = -1;
                  continue;
                }

              const unsigned int component_root =
                connected_faces.find_root(face_index);
              const auto block_id =
                component_root_to_block_id.find(component_root);

              if (block_id != component_root_to_block_id.end())
                all_faces[face_index].block_id = static_cast<int>(block_id->second);
              else
                all_faces[face_index].block_id = -1;
            }

          // 11. Measure area fragmentation as F_block = 1 - A_largest / A_total.
          // It is zero when all continental area belongs to one retained block
          // and increases as area is distributed among separate blocks.
          if (global_cont_area > 0.0)
            {
              fragmentation_index =
                1.0 - global_largest_block / global_cont_area;
            }

          if constexpr (dim == 3)
            {
              // 12. In 3D, reconstruct the coastline from mesh edges separating
              // continental and non-continental faces. Sum their chord lengths
              // to estimate the total continental perimeter.
              std::map<EdgeKey, EdgeRecord> surface_edges;

              for (unsigned int face_index = 0;
                   face_index < number_of_faces;
                   ++face_index)
                {
                  const std::vector<Point<dim>> &face_vertices =
                    all_faces[face_index].vertices;

                  for (unsigned int edge_index = 0;
                       edge_index < GeometryInfo<dim-1>::lines_per_cell;
                       ++edge_index)
                    {
                      const unsigned int first_vertex_index =
                        GeometryInfo<dim-1>::line_to_cell_vertices(edge_index, 0);
                      const unsigned int second_vertex_index =
                        GeometryInfo<dim-1>::line_to_cell_vertices(edge_index, 1);

                      std::string first_vertex_key =
                        point_key(face_vertices[first_vertex_index]);
                      std::string second_vertex_key =
                        point_key(face_vertices[second_vertex_index]);
                      if (second_vertex_key < first_vertex_key)
                        std::swap(first_vertex_key, second_vertex_key);

                      EdgeRecord &surface_edge =
                        surface_edges[ {first_vertex_key, second_vertex_key}];
                      surface_edge.incident_face_indices.push_back(face_index);
                      surface_edge.length =
                        point_distance(face_vertices[first_vertex_index],
                                       face_vertices[second_vertex_index]);
                    }
                }

              for (const auto &edge_entry : surface_edges)
                {
                  const EdgeRecord &surface_edge = edge_entry.second;

                  bool touches_continent = false;
                  bool touches_non_continent = false;

                  for (const unsigned int face_index : surface_edge.incident_face_indices)
                    {
                      touches_continent =
                        touches_continent || all_faces[face_index].is_continent;
                      touches_non_continent =
                        touches_non_continent || !all_faces[face_index].is_continent;
                    }

                  if (touches_continent && touches_non_continent)
                    global_perimeter += surface_edge.length;
                }
            }

          if (global_cont_area > 0.0)
            {
              // 13. Calculate perimeter-based fragmentation measures. P/A is
              // boundary length per unit continental area and increases for
              // smaller, narrower, or more irregular continents.
              perimeter_to_area_ratio = global_perimeter / global_cont_area;

              // Compare the measured perimeter with a planar equal-area circle:
              // F_perim = max(0, 1 - 2 sqrt(pi A) / P). Larger values represent
              // less compact outlines. On a spherical surface this reference
              // perimeter is an approximation.
              if (global_perimeter > 0.0)
                {
                  const double minimum_circle_perimeter =
                    2.0 * std::sqrt(numbers::PI * global_cont_area);
                  normalized_perimeter_fragmentation =
                    std::max(0.0,
                             1.0 - minimum_circle_perimeter / global_perimeter);
                }
            }

          std::vector<BlockSummary> blocks(global_n_blocks);
          std::vector<Point<dim>> block_centroids(global_n_blocks);

          for (unsigned int block_index = 0;
               block_index < global_n_blocks;
               ++block_index)
            blocks[block_index].id = block_index;

          // 14. Characterize each retained block by area, area-weighted drift
          // speed, and area-weighted Cartesian centroid.
          for (const FaceRecord &surface_face : all_faces)
            {
              if (surface_face.block_id < 0)
                continue;

              BlockSummary &block =
                blocks[static_cast<unsigned int>(surface_face.block_id)];
              block.area += surface_face.area;
              block.speed_area_integral += surface_face.speed * surface_face.area;
              block.n_faces += 1;

              for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
                block.centroid_numerator[coordinate] +=
                  surface_face.center[coordinate] * surface_face.area;
            }

          for (unsigned int block_index = 0;
               block_index < global_n_blocks;
               ++block_index)
            {
              if (blocks[block_index].area > 0.0)
                {
                  for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
                    block_centroids[block_index][coordinate] =
                      blocks[block_index].centroid_numerator[coordinate]
                      / blocks[block_index].area;
                }
            }

          // 15. Write optional surface and block diagnostics. The scalar results
          // are subsequently published to ASPECT's statistics and screen output.
          if (write_surface_map)
            write_surface_file(all_faces);

          if (write_block_summary)
            write_block_file(blocks, block_centroids);
        }

      // Share the global diagnostics with all MPI processes before adding them
      // to ASPECT's statistics table.
      global_cont_area =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_cont_area, 0);
      global_cont_speed_integral =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_cont_speed_integral, 0);
      global_perimeter =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_perimeter, 0);
      global_largest_block =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_largest_block, 0);
      global_n_blocks =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_n_blocks, 0);
      fragmentation_index =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), fragmentation_index, 0);
      perimeter_to_area_ratio =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), perimeter_to_area_ratio, 0);
      normalized_perimeter_fragmentation =
        Utilities::MPI::broadcast(this->get_mpi_communicator(),
                                  normalized_perimeter_fragmentation,
                                  0);

      const double mean_cont_speed_m_per_s =
        (global_cont_area > 0.0 ? global_cont_speed_integral / global_cont_area : 0.0);

      const double mean_cont_speed =
        mean_cont_speed_m_per_s * year_in_seconds;

      // Publish the global geometric and kinematic diagnostics.
      statistics.add_value("Continental area", global_cont_area);
      statistics.set_precision("Continental area", 8);
      statistics.set_scientific("Continental area", true);

      statistics.add_value("Continental perimeter", global_perimeter);
      statistics.set_precision("Continental perimeter", 8);
      statistics.set_scientific("Continental perimeter", true);

      statistics.add_value("Number of continental blocks", global_n_blocks);
      statistics.set_precision("Number of continental blocks", 0);

      statistics.add_value("Largest continental block area", global_largest_block);
      statistics.set_precision("Largest continental block area", 8);
      statistics.set_scientific("Largest continental block area", true);

      statistics.add_value("Largest-block fragmentation index", fragmentation_index);
      statistics.set_precision("Largest-block fragmentation index", 6);

      statistics.add_value("Continental perimeter-to-area ratio", perimeter_to_area_ratio);
      statistics.set_precision("Continental perimeter-to-area ratio", 8);
      statistics.set_scientific("Continental perimeter-to-area ratio", true);

      statistics.add_value("Normalized perimeter fragmentation index",
                           normalized_perimeter_fragmentation);
      statistics.set_precision("Normalized perimeter fragmentation index", 6);

      statistics.add_value("Average continental drift speed", mean_cont_speed);
      statistics.set_precision("Average continental drift speed", 8);
      statistics.set_scientific("Average continental drift speed", true);

      std::ostringstream out;
      if (output_verbose_screen_line)
        {
          out << "A_cont=" << global_cont_area
              << " m^2, P_cont=" << global_perimeter
              << " m, N_blocks=" << global_n_blocks
              << ", F_block=" << fragmentation_index
              << ", P/A=" << perimeter_to_area_ratio << " 1/m"
              << ", F_perim=" << normalized_perimeter_fragmentation
              << ", U_cont=" << mean_cont_speed << " m/yr";
        }

      return std::make_pair("Continental fragmentation statistics", out.str());
    }



    ASPECT_REGISTER_POSTPROCESSOR(
      ContinentalFragmentationStatistics,
      "continental fragmentation statistics",
      "Computes continent-specific diagnostics on the top surface from a "
      "user-defined set of compositional fields. Outputs continental area, "
      "continent perimeter, number of connected continental blocks, largest "
      "continental block area, largest-block fragmentation index, continental "
      "perimeter-to-area ratio, normalized perimeter fragmentation index, and "
      "average continental drift speed. Can also write per-face and per-block "
      "debug files.")
  }
}
