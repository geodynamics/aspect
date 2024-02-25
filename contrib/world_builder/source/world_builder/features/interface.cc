/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

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

#include "world_builder/features/interface.h"

#include <algorithm>

#include "world_builder/types/array.h"
#include "world_builder/types/object.h"
#include "world_builder/types/point.h"
#include "world_builder/world.h"

namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace Internal
    {



      /**
       * A function to turn strings into interpolation type enums.
       */
      static InterpolationType string_to_interpolation_type (const std::string &string)
      {
        if (string == "continuous monotone spline")
          {
            return InterpolationType::ContinuousMonotoneSpline;
          }

        WBAssertThrow(false,
                      "You provided an interpolation type which is not supported: " << string
                      << "This may be due to all options besides continuous monotone spline have been "
                      << "removed since version 0.5. It is best to remove the interpolation variable "
                      << "from you input file as it may be removed in future versions.");

        return InterpolationType::Invalid;
      }
    } // namespace Internal

    Interface::Interface()
      = default;

    Interface::~Interface ()
      = default;

    void
    Interface::declare_entries(Parameters &prm, const std::string &parent_name, const std::vector<std::string> &required_entries)
    {

      {
        using namespace rapidjson;
        Document &declarations = prm.declarations;

        prm.enter_subsection("defaultSnippets");
        const std::string path = prm.get_full_json_path();

        unsigned int idx = 0;
        for (auto  &it : get_snippet_map())
          {
            std::string item =  path + "/"+std::to_string(idx);

            Pointer((item).c_str()).Set(declarations, "object");
            Pointer((item+"/label").c_str()).Set(declarations,("add a '" + it.first + "'").c_str());
            prm.enter_subsection(std::to_string(idx).c_str());
            it.second(prm);
            prm.leave_subsection();

            ++idx;
          }

        prm.leave_subsection();
      }

      unsigned int counter = 0;
      for (auto  &it : get_declare_map())
        {
          prm.enter_subsection("oneOf");
          {
            prm.enter_subsection(std::to_string(counter));
            {

              prm.enter_subsection("properties");
              {
                prm.declare_entry("", Types::Object(required_entries), "feature object");

                prm.declare_entry("model", Types::String("",it.first),
                                  "The model name of the feature determining its type.");
                prm.declare_entry("name", Types::String(""),
                                  "The name which the user has given to the feature. "
                                  "This is mostly used for documentation purposes, and should in most cases be unique, "
                                  "although this is not enforced.");
                prm.declare_entry("tag", Types::String(""),
                                  "A tag which can be given to a feature. This is meant to categorize different features. "
                                  "If the tag is not provided or empty, it is set to the model name.");
                prm.declare_entry("coordinates", Types::Array(Types::Point<2>(), 1),
                                  "An array of 2d Points representing an array of coordinates where the feature is located.");

                prm.declare_entry("interpolation",Types::String("global"),
                                  "What type of interpolation should be used to enforce the minimum points per "
                                  "distance parameter. Options are 'global' and "
                                  "'continuous monotone spline' interpolation. If this "
                                  "value is set to global, the global value for interpolation is used. "
                                  "This option is deprecated and will be removed in a future release.");
                WBAssert(it.second != NULL, "No declare entries given.");
                it.second(prm, parent_name, {});
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();

          counter++;

        }
    }


    void
    Interface::get_coordinates(const std::string & /*unused*/,
                               Parameters &prm,
                               const CoordinateSystem coordinate_system)
    {
      coordinates = prm.get_vector<Point<2> >("coordinates");
      if (coordinate_system == CoordinateSystem::spherical)
        std::transform(coordinates.begin(),coordinates.end(), coordinates.begin(),
                       [](const WorldBuilder::Point<2> &p) -> WorldBuilder::Point<2> { return p *Consts::PI / 180.0;});


      // If global is given, we use the global interpolation setting, otherwise use the provided value.
      const std::string interpolation_type_string = prm.get<std::string>("interpolation") == "global" ? this->world->interpolation : prm.get<std::string>("interpolation");
      interpolation_type = WorldBuilder::Features::Internal::string_to_interpolation_type(interpolation_type_string);

      original_number_of_coordinates = coordinates.size();

      WBAssert(interpolation_type == WorldBuilder::Utilities::InterpolationType::Linear ||
               interpolation_type == WorldBuilder::Utilities::InterpolationType::MonotoneSpline ||
               interpolation_type == WorldBuilder::Utilities::InterpolationType::ContinuousMonotoneSpline,
               "For interpolation, linear and monotone spline are the only allowed values. "
               << "You provided " << interpolation_type_string << '.');

      bezier_curve = Objects::BezierCurve(coordinates);

    }


    void
    Interface::registerType(const std::string &name,
                            void ( *declare_entries)(Parameters &, const std::string &,const std::vector<std::string> &),
                            void ( *make_snippet)(Parameters &),
                            ObjectFactory *factory)
    {
      get_factory_map()[name] = factory;
      get_declare_map()[name] = declare_entries;
      get_snippet_map()[name] = make_snippet;
    }

    std::unique_ptr<Interface>
    Interface::create(const std::string &name, WorldBuilder::World *world)
    {
      std::string lower_case_name;
      std::transform(name.begin(),
                     name.end(),
                     std::back_inserter(lower_case_name),
                     ::tolower);;

      // Have a nice assert message to check whether a plugin exists in the case
      // of a debug compilation.
      WBAssertThrow(get_factory_map().find(lower_case_name) != get_factory_map().end(),
                    "Internal error: Plugin with name '" << lower_case_name << "' is not found. "
                    "The size of factories is " << get_factory_map().size() << '.');

      // Using at() because the [] will just insert values
      // which is undesirable in this case. An exception is
      // thrown when the name is not present.
      return get_factory_map().at(lower_case_name)->create(world);
    }

    Objects::PlaneDistances
    Interface::distance_to_feature_plane(const Point<3> & /*unused*/,
                                         const Objects::NaturalCoordinate & /*unused*/,
                                         const double /*unused*/) const
    {
      WBAssertThrow(false, "The distance_to_feature_plane is not yet implemented for the desinated object");
    }

  } // namespace Features
} // namespace WorldBuilder

