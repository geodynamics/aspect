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

#include <aspect/world_builder/coordinate_system/cartesian.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string.hpp>

using dealii::StandardExceptions::ExcMessage;

namespace aspect
{
  namespace WorldBuilder
  {
    namespace CoordinateSystem
    {
      Cartesian::Cartesian()
      {}

      // todo: add relative path somehow, to output when there are errors
      void
      Cartesian::read(ptree &/*tree*/)
      {
        /*boost::optional<std::string> value  = tree.get_optional<std::string> ("name");
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
            coordinates.push_back(tmp);
          }
        AssertThrow (coordinates.size() > 2, ExcMessage("This feature requires at least 3 coordinates, but only " +
                                                        dealii::Utilities::to_string(coordinates.size()) +
                                                        " where provided."));

        // Temperature submodule parameters
        value  = tree.get_optional<std::string> ("temperature submodule.name");
        AssertThrow (value, ExcMessage("Entry undeclared:  temperature submodule.name"));
        temperature_submodule_name = boost::algorithm::to_lower_copy(value.get());
        boost::algorithm::trim(temperature_submodule_name);

        if (composition_submodule_name == "constant")
          {
            value  = tree.get_optional<std::string> ("temperature submodule.depth");
            AssertThrow (value, ExcMessage("Entry undeclared:  temperature submodule.depth"));
            temperature_submodule_depth = value.get();

            value  = tree.get_optional<std::string> ("temperature submodule.temperature");
            AssertThrow (value, ExcMessage("Entry undeclared:  temperature submodule.temperature"));
            temperature_submodule_temperature = value.get();
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
            composition_submodule_depth = value.get();

            value  = tree.get_optional<std::string> ("composition submodule.composition");
            AssertThrow (value, ExcMessage("Entry undeclared:  composition submodule.temperature"));
            composition_submodule_temperature = value.get();
          }*/


      }


      aspect::WorldBuilder::Utilities::Coordinates::CoordinateSystem
      Cartesian::natural_coordinate_system() const
      {
        return aspect::WorldBuilder::Utilities::Coordinates::CoordinateSystem::cartesian;
      }


      std::array<double,3>
      Cartesian::cartesian_to_natural_coordinates(const std::array<double,3> &position) const
      {
        return position;
      }


      std::array<double,3>
      Cartesian::natural_to_cartesian_coordinates(const std::array<double,3> &position) const
      {
        return position;
      }
    }
  }
}
