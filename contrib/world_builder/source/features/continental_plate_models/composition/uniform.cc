/*
  Copyright (C) 2018 - 2021 by the authors of the World Builder code.

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

#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/nan.h>
#include <world_builder/parameters.h>

#include <world_builder/types/array.h>
#include <world_builder/types/double.h>
#include <world_builder/types/string.h>
#include <world_builder/types/object.h>
#include <world_builder/types/unsigned_int.h>
#include <world_builder/features/continental_plate_models/composition/uniform.h>


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace ContinentalPlateModels
    {
      namespace Composition
      {
        Uniform::Uniform(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          operation("")
        {
          this->world = world_;
          this->name = "uniform";
        }

        Uniform::~Uniform()
          = default;

        void
        Uniform::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Add compositions to the required parameters.
          prm.declare_entry("", Types::Object({"compositions"}), "Uniform compositional model object");


          prm.declare_entry("min depth", Types::Double(0),
                            "The depth in meters from which the composition of this feature is present.");
          prm.declare_entry("max depth", Types::Double(std::numeric_limits<double>::max()),
                            "The depth in meters to which the composition of this feature is present.");
          prm.declare_entry("compositions", Types::Array(Types::UnsignedInt(),0),
                            "A list with the labels of the composition which are present there.");
          prm.declare_entry("fractions", Types::Array(Types::Double(1.0),1),
                            "TA list of compositional fractions corresponding to the compositions list.");
          prm.declare_entry("operation", Types::String("replace",std::vector<std::string> {"replace"}),
                            "Whether the value should replace any value previously defined at this location (replace) or "
                            "add the value to the previously define value (add, not implemented). Replacing implies that all values not "
                            "explicitly defined are set to zero.");

        }

        void
        Uniform::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min depth");
          max_depth = prm.get<double>("max depth");
          compositions = prm.get_vector<unsigned int>("compositions");
          fractions = prm.get_vector<double>("fractions");
          operation = prm.get<std::string>("operation");

          WBAssertThrow(compositions.size() == fractions.size(),
                        "There are not the same amount of compositions and fractions.");
        }


        double
        Uniform::get_composition(const Point<3> & /*position*/,
                                 const double depth,
                                 const unsigned int composition_number,
                                 double composition_,
                                 const double  /*feature_min_depth*/,
                                 const double  /*feature_max_depth*/) const
        {
          double composition = composition_;
          if (depth <= max_depth && depth >= min_depth)
            {
              for (unsigned int i =0; i < compositions.size(); ++i)
                {
                  if (compositions[i] == composition_number)
                    {
                      return fractions[i];
                    }
                }

              if (operation == "replace")
                return 0.0;
            }
          return composition;
        }
        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_COMPOSITION_MODEL(Uniform, uniform)
      }
    }
  }
}
