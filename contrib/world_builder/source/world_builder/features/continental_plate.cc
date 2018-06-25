/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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

#include <aspect/world_builder/features/continental_plate.h>
#include <aspect/world_builder/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string.hpp>


using dealii::StandardExceptions::ExcMessage;

namespace aspect
{
  namespace WorldBuilder
  {
    namespace Features
    {
      ContinentalPlate::ContinentalPlate(WorldBuilder::World *world_)
        :
        temperature_submodule_depth(std::numeric_limits<double>::signaling_NaN()),
        temperature_submodule_temperature(std::numeric_limits<double>::signaling_NaN()),
        composition_submodule_depth(std::numeric_limits<double>::signaling_NaN()),
        composition_submodule_composition(std::numeric_limits<unsigned int>::signaling_NaN())
      {
        this->world = world_;
      }

      // todo: add relative path somehow, to output when there are erros
      void
      ContinentalPlate::read(ptree &tree)
      {
        boost::optional<std::string> value  = tree.get_optional<std::string> ("name");
        AssertThrow (value, ExcMessage("Entry undeclared:  name"));
        name = boost::algorithm::to_lower_copy(value.get());
        boost::algorithm::trim(name);


        boost::optional<ptree &> child = tree.get_child("coordinates");
        AssertThrow (child, ExcMessage("Entry undeclared:  coordinates"));
        for (boost::property_tree::ptree::iterator it = child.get().begin(); it != child.get().end(); ++it)
          {
            std::vector<double> tmp;
            boost::optional<ptree &> child2 = it->second.get_child("");
            AssertThrow (child, ExcMessage("This should be a 2d array, but only one dimension found."));
            for (boost::property_tree::ptree::iterator it2 = child2.get().begin(); it2 != child2.get().end(); ++it2)
              {
                tmp.push_back(dealii::Utilities::string_to_double(it2->second.get<std::string>("")));
              }
            AssertThrow (tmp.size() == 2, ExcMessage("These represent 2d coordinates, but there are " +
                                                     dealii::Utilities::to_string(tmp.size()) +
                                                     " coordinates specified."));

            std::array<double,2> tmp_array;
            std::copy(tmp.begin(), tmp.end(), tmp_array.begin());
            coordinates.push_back(tmp_array);
          }
        AssertThrow (coordinates.size() > 2, ExcMessage("This feature requires at least 3 coordinates, but only " +
                                                        dealii::Utilities::to_string(coordinates.size()) +
                                                        " where provided."));

        // Temperature submodule parameters
        value  = tree.get_optional<std::string> ("temperature submodule.name");
        AssertThrow (value, ExcMessage("Entry undeclared:  temperature submodule.name"));
        temperature_submodule_name = boost::algorithm::to_lower_copy(value.get());
        boost::algorithm::trim(temperature_submodule_name);


        if (temperature_submodule_name == "constant")
          {
            value  = tree.get_optional<std::string> ("temperature submodule.depth");
            AssertThrow (value, ExcMessage("Entry undeclared:  temperature submodule.depth"));
            temperature_submodule_depth = dealii::Utilities::string_to_double(value.get());

            value  = tree.get_optional<std::string> ("temperature submodule.temperature");
            AssertThrow (value, ExcMessage("Entry undeclared:  temperature submodule.temperature"));
            temperature_submodule_temperature = dealii::Utilities::string_to_double(value.get());
          }

        //Composition submodule parameters
        value  = tree.get_optional<std::string> ("composition submodule.name");
        AssertThrow (value, ExcMessage("Entry undeclared:  composition submodule.name"));
        composition_submodule_name = boost::algorithm::to_lower_copy(value.get());
        boost::algorithm::trim(composition_submodule_name);

        if (composition_submodule_name == "constant")
          {
            value  = tree.get_optional<std::string> ("composition submodule.depth");
            AssertThrow (value, ExcMessage("Entry undeclared:  composition submodule.depth"));
            composition_submodule_depth = dealii::Utilities::string_to_double(value.get());

            value  = tree.get_optional<std::string> ("composition submodule.composition");
            AssertThrow (value, ExcMessage("Entry undeclared:  composition submodule.temperature"));
            composition_submodule_composition = dealii::Utilities::string_to_int(value.get());
          }


      }

      double
      ContinentalPlate::temperature(const std::array<double,3> position,
                                    const double depth,
                                    const double /*gravity*/,
                                    double temperature) const
      {
        if (temperature_submodule_name == "constant")
          {
            WorldBuilder::Utilities::NaturalCoordinate natural_coordinate = WorldBuilder::Utilities::NaturalCoordinate(position,*(world->get_coordinate_system()));
            // The constant temperature module should be used for this.
            if (depth <= temperature_submodule_depth &&
                Utilities::polygon_contains_point(coordinates, natural_coordinate.get_surface_coordinates()))
              {
                // We are in the the area where the contintal plate is defined. Set the constant temperature.
                return temperature_submodule_temperature;
              }

          }
        else if (temperature_submodule_name == "none")
          {
        	return temperature;
          }
        else
          {
            AssertThrow(false,ExcMessage("Given temperature module does not exist: " + temperature_submodule_name))
          }

        return temperature;
      }

      double
      ContinentalPlate::composition(const std::array<double,3> position,
                                    const double depth,
                                    const unsigned int composition_number,
                                    double composition) const
      {
        if (composition_submodule_name == "constant")
          {
            WorldBuilder::Utilities::NaturalCoordinate natural_coordinate = WorldBuilder::Utilities::NaturalCoordinate(position,*(world->get_coordinate_system()));
            // The constant temperature module should be used for this.
            if (depth <= composition_submodule_depth &&
                Utilities::polygon_contains_point(coordinates, natural_coordinate.get_surface_coordinates()))
              {
                // We are in the the area where the contintal plate is defined. Set the constant temperature.
            	if(composition_submodule_composition == composition_number)
            	{
                return 1.0;
            	}
              }

          }
        else if (composition_submodule_name == "none")
          {
        	return composition;
          }
        else
          {
            AssertThrow(false,ExcMessage("Given composition module does not exist: " + composition_submodule_name))
          }

        return composition;
      }
    }
  }
}
