/*
  Copyright (C) 2018 - 2020 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <world_builder/features/subducting_plate.h>
#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/nan.h>
#include <world_builder/parameters.h>

#include <world_builder/types/array.h>
#include <world_builder/types/double.h>
#include <world_builder/types/point.h>
#include <world_builder/types/string.h>
#include <world_builder/types/object.h>
#include <world_builder/types/plugin_system.h>
#include <world_builder/types/unsigned_int.h>

#include <rapidjson/istreamwrapper.h>
#include "rapidjson/pointer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/error/en.h"

#include "glm/glm.h"



namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    SubductingPlate::SubductingPlate(WorldBuilder::World *world_)
      :
      reference_point(0,0,cartesian)
    {
      this->world = world_;
      this->name = "subducting plate";
    }

    SubductingPlate::~SubductingPlate()
    { }



    void
    SubductingPlate::declare_entries(Parameters &prm,
                                     const std::string &parent_name,
                                     const std::vector<std::string> & /*required_entries*/)
    {
      // This statment is needed because of the recursion associated with
      // the sections entry.
      if (parent_name == "items")
        prm.enter_subsection("properties");

      prm.declare_entry("min depth", Types::Double(0),
                        "The depth to which this feature is present");
      prm.declare_entry("max depth", Types::Double(std::numeric_limits<double>::max()),
                        "The depth to which this feature is present");
      prm.declare_entry("dip point", Types::Point<2>(),
                        "The depth to which this feature is present");

      /*prm.declare_entry("segments", Types::Array(Types::Segment(0,Point<2>(0,0,invalid),Point<2>(0,0,invalid),Point<2>(0,0,invalid),
                                                                Types::PluginSystem("", Features::SubductingPlateModels::Temperature::Interface::declare_entries, {"model"}),
                                                                Types::PluginSystem("", Features::SubductingPlateModels::Composition::Interface::declare_entries, {"model"}),
                                                                Types::PluginSystem("", Features::SubductingPlateModels::Grains::Interface::declare_entries, {"model"}))),
                        "The depth to which this feature is present");*/
      prm.declare_entry("segments", Types::Array(Types::Segment(0,Point<2>(0,0,invalid),Point<2>(0,0,invalid),Point<2>(0,0,invalid),
                                                                Types::PluginSystem("", Features::SubductingPlateModels::Temperature::Interface::declare_entries, {"model"}),
                                                                Types::PluginSystem("", Features::SubductingPlateModels::Composition::Interface::declare_entries, {"model"}),
                                                                Types::PluginSystem("", Features::SubductingPlateModels::Grains::Interface::declare_entries, {"model"}))),
                        "The depth to which this feature is present");

      prm.declare_entry("temperature models",
                        Types::PluginSystem("", Features::SubductingPlateModels::Temperature::Interface::declare_entries, {"model"}),
                        "A list of temperature models.");
      prm.declare_entry("composition models",
                        Types::PluginSystem("", Features::SubductingPlateModels::Composition::Interface::declare_entries, {"model"}),
                        "A list of composition models.");
      prm.declare_entry("grains models",
                        Types::PluginSystem("", Features::SubductingPlateModels::Grains::Interface::declare_entries, {"model"}),
                        "A list of grains models.");

      if (parent_name != "items")
        {
          // This only happens if we are not in sections
          prm.declare_entry("sections", Types::Array(Types::PluginSystem("",Features::SubductingPlate::declare_entries, {"coordinate"}, false)),"A list of feature properties for a coordinate.");
        }
      else
        {

          // this only happens in sections
          prm.declare_entry("coordinate", Types::UnsignedInt(0),
                            "The coordinate which should be overwritten");

          prm.leave_subsection();
        }
    }

    void
    SubductingPlate::parse_entries(Parameters &prm)
    {
      const CoordinateSystem coordinate_system = prm.coordinate_system->natural_coordinate_system();

      this->name = prm.get<std::string>("name");
      this->get_coordinates("coordinates", prm, coordinate_system);


      starting_depth = prm.get<double>("min depth");
      maximum_depth = prm.get<double>("max depth");

      const size_t n_sections = this->original_number_of_coordinates;

      reference_point = prm.get<Point<2> >("dip point");

      default_temperature_models.resize(0);
      default_composition_models.resize(0);
      default_grains_models.resize(0);
      prm.get_shared_pointers<Features::SubductingPlateModels::Temperature::Interface>("temperature models", default_temperature_models);
      prm.get_shared_pointers<Features::SubductingPlateModels::Composition::Interface>("composition models", default_composition_models);
      prm.get_shared_pointers<Features::SubductingPlateModels::Grains::Interface>("grains models", default_grains_models);

      // get the default segments.
      default_segment_vector = prm.get_vector<Objects::Segment<Features::SubductingPlateModels::Temperature::Interface,
      Features::SubductingPlateModels::Composition::Interface,
      Features::SubductingPlateModels::Grains::Interface> >("segments", default_temperature_models, default_composition_models, default_grains_models);


      // This vector stores segments to this coordiante/section.
      //First used (raw) pointers to the segment relevant to this coordinate/section,
      // but I do not trust it won't fail when memory is moved. So storing the all the data now.
      segment_vector.resize(0);
      segment_vector.resize(n_sections, default_segment_vector);


      // now search whether a section is present, if so, replace the default segments.
      std::vector<std::unique_ptr<Features::SubductingPlate> > sections_vector;
      prm.get_unique_pointers("sections", sections_vector);

      prm.enter_subsection("sections");
      for (unsigned int i_section = 0; i_section < n_sections; ++i_section)
        {
          // first check whether this section/coordinate has a a special overwrite
          for (unsigned int i_sector = 0; i_sector < sections_vector.size(); ++i_sector)
            {
              prm.enter_subsection(std::to_string(i_sector));
              {
                const unsigned int change_coord_number = prm.get<unsigned int>("coordinate");

                WBAssertThrow(segment_vector.size() > change_coord_number, "Error: for subducting plate with name: '" << this->name
                              << "', trying to change the section of coordinate " << change_coord_number
                              << " while only " << segment_vector.size() << " coordinates are defined.");

                std::vector<std::shared_ptr<Features::SubductingPlateModels::Temperature::Interface> > local_default_temperature_models;
                std::vector<std::shared_ptr<Features::SubductingPlateModels::Composition::Interface>  > local_default_composition_models;
                std::vector<std::shared_ptr<Features::SubductingPlateModels::Grains::Interface>  > local_default_grains_models;

                if (prm.get_shared_pointers<Features::SubductingPlateModels::Temperature::Interface>("temperature models", local_default_temperature_models) == false)
                  {
                    // no local temperature model, use global default
                    local_default_temperature_models = default_temperature_models;
                  }

                if (prm.get_shared_pointers<Features::SubductingPlateModels::Composition::Interface>("composition models", local_default_composition_models) == false)
                  {
                    // no local composition model, use global default
                    local_default_composition_models = default_composition_models;
                  }

                if (prm.get_shared_pointers<Features::SubductingPlateModels::Grains::Interface>("grains models", local_default_grains_models) == false)
                  {
                    // no local composition model, use global default
                    local_default_grains_models = default_grains_models;
                  }

                segment_vector[change_coord_number] = prm.get_vector<Objects::Segment<Features::SubductingPlateModels::Temperature::Interface,
                                                      Features::SubductingPlateModels::Composition::Interface,
                                                      Features::SubductingPlateModels::Grains::Interface> >("segments", local_default_temperature_models, local_default_composition_models, local_default_grains_models);


                WBAssertThrow(segment_vector[change_coord_number].size() == default_segment_vector.size(),
                              "Error: There are not the same amount of segments in section with coordinate " << change_coord_number
                              << " (" << segment_vector[change_coord_number].size() << " segments) as in the default segment ("
                              << default_segment_vector.size() << " segments). This is not allowed.");

                prm.enter_subsection("segments");
                {
                  for (unsigned int i = 0; i < segment_vector[change_coord_number].size(); ++i)
                    {
                      prm.enter_subsection(std::to_string(i));
                      {
                        prm.enter_subsection("temperature models");
                        {
                          for (unsigned int j = 0; j < segment_vector[change_coord_number][i].temperature_systems.size(); ++j)
                            {
                              prm.enter_subsection(std::to_string(j));
                              {
                                segment_vector[change_coord_number][i].temperature_systems[j]->parse_entries(prm);
                              }
                              prm.leave_subsection();
                            }
                        }
                        prm.leave_subsection();


                        prm.enter_subsection("composition models");
                        {
                          for (unsigned int j = 0; j < segment_vector[change_coord_number][i].composition_systems.size(); ++j)
                            {
                              prm.enter_subsection(std::to_string(j));
                              {
                                segment_vector[change_coord_number][i].composition_systems[j]->parse_entries(prm);
                              }
                              prm.leave_subsection();
                            }
                        }
                        prm.leave_subsection();

                        prm.enter_subsection("grains models");
                        {
                          for (unsigned int j = 0; j < segment_vector[change_coord_number][i].grains_systems.size(); ++j)
                            {
                              prm.enter_subsection(std::to_string(j));
                              {
                                segment_vector[change_coord_number][i].grains_systems[j]->parse_entries(prm);
                              }
                              prm.leave_subsection();
                            }
                        }
                        prm.leave_subsection();
                      }
                      prm.leave_subsection();
                    }
                }
                prm.leave_subsection();

              }
              prm.leave_subsection();
            }

        }
      prm.leave_subsection();


      prm.enter_subsection("segments");
      {
        for (unsigned int i = 0; i < default_segment_vector.size(); ++i)
          {
            prm.enter_subsection(std::to_string(i));
            {
              prm.enter_subsection("temperature models");
              {
                for (unsigned int j = 0; j < default_segment_vector[i].temperature_systems.size(); ++j)
                  {
                    prm.enter_subsection(std::to_string(j));
                    {
                      default_segment_vector[i].temperature_systems[j]->parse_entries(prm);
                    }
                    prm.leave_subsection();
                  }
              }
              prm.leave_subsection();


              prm.enter_subsection("composition models");
              {
                for (unsigned int j = 0; j < default_segment_vector[i].composition_systems.size(); ++j)
                  {
                    prm.enter_subsection(std::to_string(j));
                    {
                      default_segment_vector[i].composition_systems[j]->parse_entries(prm);
                    }
                    prm.leave_subsection();
                  }
              }
              prm.leave_subsection();


              prm.enter_subsection("grains models");
              {
                for (unsigned int j = 0; j < default_segment_vector[i].grains_systems.size(); ++j)
                  {
                    prm.enter_subsection(std::to_string(j));
                    {
                      default_segment_vector[i].grains_systems[j]->parse_entries(prm);
                    }
                    prm.leave_subsection();
                  }
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();

      maximum_slab_thickness = 0;
      maximum_total_slab_length = 0;
      total_slab_length.resize(original_number_of_coordinates);
      slab_segment_lengths.resize(original_number_of_coordinates);
      slab_segment_thickness.resize(original_number_of_coordinates);
      slab_segment_top_truncation.resize(original_number_of_coordinates);
      slab_segment_angles.resize(original_number_of_coordinates);

      for (unsigned int i = 0; i < segment_vector.size(); ++i)
        {
          double local_total_slab_length = 0;
          slab_segment_lengths[i].resize(segment_vector[i].size());
          slab_segment_thickness[i].resize(segment_vector[i].size(), Point<2>(invalid));
          slab_segment_top_truncation[i].resize(segment_vector[i].size(), Point<2>(invalid));
          slab_segment_angles[i].resize(segment_vector[i].size(), Point<2>(invalid));
          for (unsigned int j = 0; j < segment_vector[i].size(); ++j)
            {
              slab_segment_lengths[i][j] = segment_vector[i][j].value_length;
              local_total_slab_length += segment_vector[i][j].value_length;

              slab_segment_thickness[i][j] = segment_vector[i][j].value_thickness;
              slab_segment_top_truncation[i][j] = segment_vector[i][j].value_top_truncation;

              slab_segment_angles[i][j] = segment_vector[i][j].value_angle * (const_pi/180);
            }
          maximum_slab_thickness = std::max(maximum_slab_thickness, local_total_slab_length);
          total_slab_length[i] = local_total_slab_length;
          maximum_total_slab_length = std::max(maximum_total_slab_length, local_total_slab_length);
        }
    }


    double
    SubductingPlate::temperature(const Point<3> &position,
                                 const double depth,
                                 const double gravity_norm,
                                 double temperature) const
    {
      WorldBuilder::Utilities::NaturalCoordinate natural_coordinate = WorldBuilder::Utilities::NaturalCoordinate(position,
                                                                      *(world->parameters.coordinate_system));

      // The depth variable is the distance from the surface to the position, the depth
      // coordinate is the distance from the bottom of the model to the position and
      // the starting radius is the distance from the bottom of the model to the surface.
      const double starting_radius = natural_coordinate.get_depth_coordinate() + depth - starting_depth;

      WBAssert(std::abs(starting_radius) > std::numeric_limits<double>::epsilon(), "World Builder error: starting_radius can not be zero. "
               << "Position = " << position[0] << ":" << position[1] << ":" << position[2]
               << ", natural_coordinate.get_depth_coordinate() = " << natural_coordinate.get_depth_coordinate()
               << ", depth = " << depth
               << ", starting_depth " << starting_depth
              );

      // todo: explain and check -starting_depth
      if (depth <= maximum_depth && depth >= starting_depth && depth <= maximum_total_slab_length + maximum_slab_thickness)
        {
          /*WBAssert(coordinates.size() == slab_segment_lengths.size(),
                   "Internal error: The size of coordinates (" << coordinates.size()
                   << ") and slab_segment_lengths (" << slab_segment_lengths.size() << ") are different.");
          WBAssert(coordinates.size() == slab_segment_angles.size(),
                   "Internal error: The size of coordinates (" << coordinates.size()
                   << ") and slab_segment_angles (" << slab_segment_angles.size() << ") are different.");
          WBAssert(coordinates.size() == slab_segment_angles.size(),
                   "Internal error: The size of coordinates (" << coordinates.size()
                   << ") and one_dimensional_coordinates (" << one_dimensional_coordinates.size() << ") are different.");*/
          // todo: explain
          std::map<std::string,double> distance_from_planes =
            WorldBuilder::Utilities::distance_point_from_curved_planes(position,
                                                                       reference_point,
                                                                       coordinates,
                                                                       slab_segment_lengths,
                                                                       slab_segment_angles,
                                                                       starting_radius,
                                                                       this->world->parameters.coordinate_system,
                                                                       false,
                                                                       one_dimensional_coordinates);

          const double distance_from_plane = distance_from_planes["distanceFromPlane"];
          const double distance_along_plane = distance_from_planes["distanceAlongPlane"];
          const double section_fraction = distance_from_planes["sectionFraction"];
          const size_t current_section = static_cast<size_t>(std::floor(one_dimensional_coordinates[static_cast<size_t>(distance_from_planes["section"])]));
          const size_t next_section = current_section + 1;
          const size_t current_segment = static_cast<size_t>(distance_from_planes["segment"]); // the original value was a unsigned in, converting it back.
          //const size_t next_segment = current_segment + 1;
          const double segment_fraction = distance_from_planes["segmentFraction"];

          if (abs(distance_from_plane) < INFINITY || (distance_along_plane) < INFINITY)
            {
              // We want to do both section (horizontal) and segment (vertical) interpolation.
              // first for thickness
              const double thickness_up = slab_segment_thickness[current_section][current_segment][0]
                                          + section_fraction
                                          * (slab_segment_thickness[next_section][current_segment][0]
                                             - slab_segment_thickness[current_section][current_segment][0]);
              const double thickness_down = slab_segment_thickness[current_section][current_segment][1]
                                            + section_fraction
                                            * (slab_segment_thickness[next_section][current_segment][1]
                                               - slab_segment_thickness[current_section][current_segment][1]);
              const double thickness_local = thickness_up + segment_fraction * (thickness_down - thickness_up);
              distance_from_planes["thicknessLocal"] = thickness_local;

              // secondly for top truncation
              const double top_truncation_up = slab_segment_top_truncation[current_section][current_segment][0]
                                               + section_fraction
                                               * (slab_segment_top_truncation[next_section][current_segment][0]
                                                  - slab_segment_top_truncation[current_section][current_segment][0]);
              const double top_truncation_down = slab_segment_top_truncation[current_section][current_segment][1]
                                                 + section_fraction
                                                 * (slab_segment_top_truncation[next_section][current_segment][1]
                                                    - slab_segment_top_truncation[current_section][current_segment][1]);
              const double top_truncation_local = top_truncation_up + segment_fraction * (top_truncation_down - top_truncation_up);

              // if the thickness is zero, we don't need to compute anything, so return.
              if (std::fabs(thickness_local) < 2.0 * std::numeric_limits<double>::epsilon())
                return temperature;

              // if the thickness is smaller than what is truncated off at the top, we don't need to compute anything, so return.
              if (thickness_local < top_truncation_local)
                return temperature;

              const double max_slab_length = total_slab_length[current_section] +
                                             section_fraction *
                                             (total_slab_length[next_section] - total_slab_length[current_section]);

              if (distance_from_plane >= top_truncation_local &&
                  distance_from_plane <= thickness_local &&
                  distance_along_plane >= 0 &&
                  distance_along_plane <= max_slab_length)
                {
                  // Inside the slab!
                  double temperature_current_section = temperature;
                  double temperature_next_section = temperature;

                  for (auto &temperature_model: segment_vector[current_section][current_segment].temperature_systems)
                    {
                      temperature_current_section = temperature_model->get_temperature(position,
                                                                                       depth,
                                                                                       gravity_norm,
                                                                                       temperature_current_section,
                                                                                       starting_depth,
                                                                                       maximum_depth,
                                                                                       distance_from_planes);

                      WBAssert(!std::isnan(temperature_current_section), "Temparture is not a number: " << temperature_current_section
                               << ", based on a temperature model with the name " << temperature_model->get_name());
                      WBAssert(std::isfinite(temperature_current_section), "Temparture is not a finite: " << temperature_current_section
                               << ", based on a temperature model with the name " << temperature_model->get_name());

                    }

                  for (auto &temperature_model: segment_vector[next_section][current_segment].temperature_systems)
                    {
                      temperature_next_section = temperature_model->get_temperature(position,
                                                                                    depth,
                                                                                    gravity_norm,
                                                                                    temperature_next_section,
                                                                                    starting_depth,
                                                                                    maximum_depth,
                                                                                    distance_from_planes);

                      WBAssert(!std::isnan(temperature_next_section), "Temparture is not a number: " << temperature_next_section
                               << ", based on a temperature model with the name " << temperature_model->get_name());
                      WBAssert(std::isfinite(temperature_next_section), "Temparture is not a finite: " << temperature_next_section
                               << ", based on a temperature model with the name " << temperature_model->get_name());

                    }

                  // linear interpolation between current and next section temperatures
                  temperature = temperature_current_section + section_fraction * (temperature_next_section - temperature_current_section);
                }
            }
        }
      return temperature;
    }

    double
    SubductingPlate::composition(const Point<3> &position,
                                 const double depth,
                                 const unsigned int composition_number,
                                 double composition) const
    {

      WorldBuilder::Utilities::NaturalCoordinate natural_coordinate = WorldBuilder::Utilities::NaturalCoordinate(position,
                                                                      *(world->parameters.coordinate_system));
      // todo: explain
      const double starting_radius = natural_coordinate.get_depth_coordinate() + depth - starting_depth;

      // todo: explain and check -starting_depth
      if (depth <= maximum_depth && depth >= starting_depth && depth <= maximum_total_slab_length + maximum_slab_thickness)
        {
          // todo: explain
          std::map<std::string,double> distance_from_planes =
            WorldBuilder::Utilities::distance_point_from_curved_planes(position,
                                                                       reference_point,
                                                                       coordinates,
                                                                       slab_segment_lengths,
                                                                       slab_segment_angles,
                                                                       starting_radius,
                                                                       this->world->parameters.coordinate_system,
                                                                       false,
                                                                       one_dimensional_coordinates);

          const double distance_from_plane = distance_from_planes["distanceFromPlane"];
          const double distance_along_plane = distance_from_planes["distanceAlongPlane"];
          const double section_fraction = distance_from_planes["sectionFraction"];
          const size_t current_section = static_cast<size_t>(std::floor(one_dimensional_coordinates[static_cast<size_t>(distance_from_planes["section"])]));
          const size_t next_section = current_section + 1;
          const size_t current_segment = static_cast<size_t>(distance_from_planes["segment"]); // the original value was a unsigned in, converting it back.
          //const size_t next_segment = current_segment + 1;
          const double segment_fraction = distance_from_planes["segmentFraction"];

          if (abs(distance_from_plane) < INFINITY || (distance_along_plane) < INFINITY)
            {
              // We want to do both section (horizontal) and segment (vertical) interpolation.

              // We want to do both section (horizontal) and segment (vertical) interpolation.
              // first for thickness
              const double thickness_up = slab_segment_thickness[current_section][current_segment][0]
                                          + section_fraction
                                          * (slab_segment_thickness[next_section][current_segment][0]
                                             - slab_segment_thickness[current_section][current_segment][0]);
              const double thickness_down = slab_segment_thickness[current_section][current_segment][1]
                                            + section_fraction
                                            * (slab_segment_thickness[next_section][current_segment][1]
                                               - slab_segment_thickness[current_section][current_segment][1]);
              const double thickness_local = thickness_up + segment_fraction * (thickness_down - thickness_up);
              distance_from_planes["thicknessLocal"] = thickness_local;

              // secondly for top truncation
              const double top_truncation_up = slab_segment_top_truncation[current_section][current_segment][0]
                                               + section_fraction
                                               * (slab_segment_top_truncation[next_section][current_segment][0]
                                                  - slab_segment_top_truncation[current_section][current_segment][0]);
              const double top_truncation_down = slab_segment_top_truncation[current_section][current_segment][1]
                                                 + section_fraction
                                                 * (slab_segment_top_truncation[next_section][current_segment][1]
                                                    - slab_segment_top_truncation[current_section][current_segment][1]);
              const double top_truncation_local = top_truncation_up + segment_fraction * (top_truncation_down - top_truncation_up);

              // if the thickness is zero, we don't need to compute anything, so return.
              if (std::fabs(thickness_local) < 2.0 * std::numeric_limits<double>::epsilon())
                return composition;

              // if the thickness is smaller than what is truncated off at the top, we don't need to compute anything, so return.
              if (thickness_local < top_truncation_local)
                return composition;

              const double max_slab_length = total_slab_length[current_section] +
                                             section_fraction *
                                             (total_slab_length[next_section] - total_slab_length[current_section]);

              if (distance_from_plane >= top_truncation_local &&
                  distance_from_plane <= thickness_local &&
                  distance_along_plane >= 0 &&
                  distance_along_plane <= max_slab_length)
                {
                  // Inside the slab!

                  double composition_current_section = composition;
                  double composition_next_section = composition;

                  for (auto &composition_model: segment_vector[current_section][current_segment].composition_systems)
                    {
                      composition_current_section = composition_model->get_composition(position,
                                                                                       depth,
                                                                                       composition_number,
                                                                                       composition_current_section,
                                                                                       starting_depth,
                                                                                       maximum_depth,
                                                                                       distance_from_planes);

                      WBAssert(!std::isnan(composition_current_section), "Composition_current_section is not a number: " << composition_current_section
                               << ", based on a temperature model with the name " << composition_model->get_name());
                      WBAssert(std::isfinite(composition_current_section), "Composition_current_section is not a finite: " << composition_current_section
                               << ", based on a temperature model with the name " << composition_model->get_name());

                    }

                  for (auto &composition_model: segment_vector[next_section][current_segment].composition_systems)
                    {
                      composition_next_section = composition_model->get_composition(position,
                                                                                    depth,
                                                                                    composition_number,
                                                                                    composition_next_section,
                                                                                    starting_depth,
                                                                                    maximum_depth,
                                                                                    distance_from_planes);

                      WBAssert(!std::isnan(composition_next_section), "Composition_next_section is not a number: " << composition_next_section
                               << ", based on a temperature model with the name " << composition_model->get_name());
                      WBAssert(std::isfinite(composition_next_section), "Composition_next_section is not a finite: " << composition_next_section
                               << ", based on a temperature model with the name " << composition_model->get_name());

                    }

                  // linear interpolation between current and next section temperatures
                  composition = composition_current_section + section_fraction * (composition_next_section - composition_current_section);


                }
            }
        }

      return composition;
    }


    WorldBuilder::grains
    SubductingPlate::grains(const Point<3> &position,
                            const double depth,
                            const unsigned int composition_number,
                            WorldBuilder::grains grains) const
    {

      WorldBuilder::Utilities::NaturalCoordinate natural_coordinate = WorldBuilder::Utilities::NaturalCoordinate(position,
                                                                      *(world->parameters.coordinate_system));
      // todo: explain
      const double starting_radius = natural_coordinate.get_depth_coordinate() + depth - starting_depth;

      // todo: explain and check -starting_depth
      if (depth <= maximum_depth && depth >= starting_depth && depth <= maximum_total_slab_length + maximum_slab_thickness)
        {
          // todo: explain
          std::map<std::string,double> distance_from_planes =
            WorldBuilder::Utilities::distance_point_from_curved_planes(position,
                                                                       reference_point,
                                                                       coordinates,
                                                                       slab_segment_lengths,
                                                                       slab_segment_angles,
                                                                       starting_radius,
                                                                       this->world->parameters.coordinate_system,
                                                                       false,
                                                                       one_dimensional_coordinates);

          const double distance_from_plane = distance_from_planes["distanceFromPlane"];
          const double distance_along_plane = distance_from_planes["distanceAlongPlane"];
          const double section_fraction = distance_from_planes["sectionFraction"];
          const size_t current_section = static_cast<size_t>(std::floor(one_dimensional_coordinates[static_cast<size_t>(distance_from_planes["section"])]));
          const size_t next_section = current_section + 1;
          const size_t current_segment = static_cast<size_t>(distance_from_planes["segment"]); // the original value was a unsigned in, converting it back.
          //const size_t next_segment = current_segment + 1;
          const double segment_fraction = distance_from_planes["segmentFraction"];

          if (abs(distance_from_plane) < INFINITY || (distance_along_plane) < INFINITY)
            {
              // We want to do both section (horizontal) and segment (vertical) interpolation.

              // We want to do both section (horizontal) and segment (vertical) interpolation.
              // first for thickness
              const double thickness_up = slab_segment_thickness[current_section][current_segment][0]
                                          + section_fraction
                                          * (slab_segment_thickness[next_section][current_segment][0]
                                             - slab_segment_thickness[current_section][current_segment][0]);
              const double thickness_down = slab_segment_thickness[current_section][current_segment][1]
                                            + section_fraction
                                            * (slab_segment_thickness[next_section][current_segment][1]
                                               - slab_segment_thickness[current_section][current_segment][1]);
              const double thickness_local = thickness_up + segment_fraction * (thickness_down - thickness_up);
              distance_from_planes["thicknessLocal"] = thickness_local;

              // secondly for top truncation
              const double top_truncation_up = slab_segment_top_truncation[current_section][current_segment][0]
                                               + section_fraction
                                               * (slab_segment_top_truncation[next_section][current_segment][0]
                                                  - slab_segment_top_truncation[current_section][current_segment][0]);
              const double top_truncation_down = slab_segment_top_truncation[current_section][current_segment][1]
                                                 + section_fraction
                                                 * (slab_segment_top_truncation[next_section][current_segment][1]
                                                    - slab_segment_top_truncation[current_section][current_segment][1]);
              const double top_truncation_local = top_truncation_up + segment_fraction * (top_truncation_down - top_truncation_up);

              // if the thickness is zero, we don't need to compute anything, so return.
              if (std::fabs(thickness_local) < 2.0 * std::numeric_limits<double>::epsilon())
                return grains;

              // if the thickness is smaller than what is truncated off at the top, we don't need to compute anything, so return.
              if (thickness_local < top_truncation_local)
                return grains;

              const double max_slab_length = total_slab_length[current_section] +
                                             section_fraction *
                                             (total_slab_length[next_section] - total_slab_length[current_section]);

              if (distance_from_plane >= top_truncation_local &&
                  distance_from_plane <= thickness_local &&
                  distance_along_plane >= 0 &&
                  distance_along_plane <= max_slab_length)
                {
                  // Inside the slab!
                  WorldBuilder::grains  grains_current_section = grains;
                  WorldBuilder::grains  grains_next_section = grains;

                  for (auto &grains_model: segment_vector[current_section][current_segment].grains_systems)
                    {
                      grains_current_section = grains_model->get_grains(position,
                                                                        depth,
                                                                        composition_number,
                                                                        grains_current_section,
                                                                        starting_depth,
                                                                        maximum_depth,
                                                                        distance_from_planes);

                      /*WBAssert(!std::isnan(composition_current_section), "Composition_current_section is not a number: " << composition_current_section
                               << ", based on a temperature model with the name " << composition_model->get_name());
                      WBAssert(std::isfinite(composition_current_section), "Composition_current_section is not a finite: " << composition_current_section
                               << ", based on a temperature model with the name " << composition_model->get_name());*/

                    }

                  for (auto &grains_model: segment_vector[next_section][current_segment].grains_systems)
                    {
                      grains_next_section = grains_model->get_grains(position,
                                                                     depth,
                                                                     composition_number,
                                                                     grains_next_section,
                                                                     starting_depth,
                                                                     maximum_depth,
                                                                     distance_from_planes);

                      /*WBAssert(!std::isnan(composition_next_section), "Composition_next_section is not a number: " << composition_next_section
                               << ", based on a temperature model with the name " << composition_model->get_name());
                      WBAssert(std::isfinite(composition_next_section), "Composition_next_section is not a finite: " << composition_next_section
                               << ", based on a temperature model with the name " << composition_model->get_name());*/

                    }

                  // linear interpolation between current and next section temperatures
                  for (size_t i = 0; i < grains.sizes.size(); i++)
                    {
                      grains.sizes[i] = grains_current_section.sizes[i] + section_fraction * (grains_next_section.sizes[i] - grains_current_section.sizes[i]);
                    }

                  // average two rotations matrices throu quaternions.
                  for (size_t i = 0; i < grains_current_section.rotation_matrices.size(); i++)
                    {
                      glm::quaternion::quat quat_current = glm::quaternion::quat_cast(grains_current_section.rotation_matrices[i]);
                      glm::quaternion::quat quat_next = glm::quaternion::quat_cast(grains_next_section.rotation_matrices[i]);

                      glm::quaternion::quat quat_average = glm::quaternion::slerp(quat_current,quat_next,section_fraction);

                      grains.rotation_matrices[i] = glm::quaternion::mat3_cast(quat_average);
                    }


                }
            }
        }

      return grains;
    }

    /**
     * Register plugin
     */
    WB_REGISTER_FEATURE(SubductingPlate, subducting plate)
  }
}

