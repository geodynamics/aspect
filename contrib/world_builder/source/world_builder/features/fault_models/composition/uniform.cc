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

#include "world_builder/features/fault_models/composition/uniform.h"


#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/unsigned_int.h"
#include "world_builder/utilities.h"


namespace WorldBuilder
{

  using namespace Utilities;

  namespace Features
  {
    using namespace FeatureUtilities;
    namespace FaultModels
    {
      namespace Composition
      {
        Uniform::Uniform(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN)
        {
          this->world = world_;
          this->name = "uniform";
        }

        Uniform::~Uniform()
          = default;

        void
        Uniform::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add compositions to the required parameters.
          prm.declare_entry("", Types::Object({"compositions"}),
                            "Uniform compositional model. Sets constant compositional field.");

          // Declare entries of this plugin
          prm.declare_entry("min distance fault center", Types::Double(0),
                            "The distance in meters from which the composition of this feature is present.");
          prm.declare_entry("max distance fault center", Types::Double(std::numeric_limits<double>::max()),
                            "The distance in meters to which the composition of this feature is present.");
          prm.declare_entry("compositions", Types::Array(Types::UnsignedInt(),0),
                            "A list with the labels of the composition which are present there.");
          prm.declare_entry("fractions", Types::Array(Types::Double(1.0),1),
                            "TA list of compositional fractions corresponding to the compositions list.");
          prm.declare_entry("operation", Types::String("replace", std::vector<std::string> {"replace", "replace defined only", "add", "subtract"}),
                            "Whether the value should replace any value previously defined at this location (replace) or "
                            "add the value to the previously define value. Replacing implies that all compositions not "
                            "explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.");

        }

        void
        Uniform::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min distance fault center");
          max_depth = prm.get<double>("max distance fault center");
          compositions = prm.get_vector<unsigned int>("compositions");
          fractions = prm.get_vector<double>("fractions");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));

          WBAssertThrow(compositions.size() == fractions.size(),
                        "There are not the same amount of compositions and fractions.");
        }


        double
        Uniform::get_composition(const Point<3> & /*position_in_cartesian_coordinates*/,
                                 const double  /*depth*/,
                                 const unsigned int composition_number,
                                 double composition,
                                 const double  /*feature_min_depth*/,
                                 const double  /*feature_max_depth*/,
                                 const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_plane,
                                 const AdditionalParameters & /*additional_parameters*/) const
        {
          if (std::fabs(distance_from_plane.distance_from_plane) <= max_depth && std::fabs(distance_from_plane.distance_from_plane) >= min_depth)
            {
              for (unsigned int i =0; i < compositions.size(); ++i)
                {
                  if (compositions[i] == composition_number)
                    {
                      return apply_operation(operation,composition,fractions[i]);
                    }
                }

              if (operation == Operations::REPLACE)
                return 0.0;
            }
          return composition;
        }
        WB_REGISTER_FEATURE_FAULT_COMPOSITION_MODEL(Uniform, uniform)
      } // namespace Composition
    } // namespace FaultModels
  } // namespace Features
} // namespace WorldBuilder


