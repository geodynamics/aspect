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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/world_builder/continental_plate.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>

namespace aspect
{
  namespace WorldBuilder
  {
    template <int dim>
    double
    ContinentalPlate<dim>::
    temperature (const Utilities::NaturalCoordinate<3> &position,
                 const unsigned int object_number,
                 const double depth,
                 double temperature) const
    {
      // inside a continental plate
      if (temperature_submodules[object_number] == "constant")
        {
          // The constant temperature module should be used for this.
          if (depth <= temperature_depths[object_number] &&
              Utilities::polygon_contains_point<3>(object_coordinates[object_number],
                                                   Utilities::convert_array_to_point<2>(
                                                     position.get_surface_coordinates())))
            {
              // We are in the the area where the contintal plate is defined. Set the constant temperature.
              return temperatures[object_number];
            }
        }
      return temperature;
    }

    template <int dim>
    double
    ContinentalPlate<dim>::
    composition (const Utilities::NaturalCoordinate<3> &position,
                 const unsigned int n_comp,
                 const unsigned int object_number,
                 const double depth,
                 double  composition) const
    {
      // inside a continental plate
      if (composition_submodules[object_number] == "constant")
        {
          // The constant temperature module should be used for this.
          if (depth <= composition_depths[object_number] &&
              Utilities::polygon_contains_point<3>(object_coordinates[object_number],
                                                   Utilities::convert_array_to_point<2>(
                                                     position.get_surface_coordinates())))
            {
              // We are in the the area where the contintal plate is defined. Set the constant composition.
              if (n_comp == compositions[object_number])
                {
                  return 1.0;
                }
              else
                {
                  return 0.0;
                }
            }
        }
      return composition;
    }

    template <int dim>
    void
    ContinentalPlate<dim>::
    declare_parameters (ParameterHandler &)
    {
      WorldBuilderParameterHandler &wb_prm = Manager<dim>::ph;
      wb_prm.enter_subsection("Surface objects");
      {
        wb_prm.enter_subsection("array");
        {
          wb_prm.enter_subsection("module");
          {
            wb_prm.declare_entry("name","", Patterns::Anything(),"TODO");
            wb_prm.enter_subsection("temperature submodule");
            {
              wb_prm.declare_entry("name","none", Patterns::Anything(),"TODO");
              wb_prm.declare_entry("temperature","273.15", Patterns::Double(),"TODO");
              wb_prm.declare_entry("depth","100e3", Patterns::Double(),"TODO");
            }
            wb_prm.leave_subsection();
            wb_prm.enter_subsection("composition submodule");
            {
              wb_prm.declare_entry("name","none", Patterns::Anything(),"TODO");
              wb_prm.declare_entry("composition","0", Patterns::Double(),"TODO");
              wb_prm.declare_entry("depth","100e3", Patterns::Anything(),"TODO");
            }
            wb_prm.leave_subsection();
          }
          wb_prm.leave_subsection();
        }
        wb_prm.leave_subsection();
      }
      wb_prm.leave_subsection();

    }

    template <int dim>
    void
    ContinentalPlate<dim>::
    parse_parameters (ParameterHandler &)
    {
      WorldBuilderParameterHandler &wb_prm = Manager<dim>::ph;
      wb_prm.enter_subsection("Surface objects");
      {
        unsigned int number_of_objects = Utilities::string_to_int(wb_prm.get("size"));
        // We allocate everything with the size of number of objects.
        // This may be a bit of a waste of memory, but I think this
        // is too little memory in total to worry about and it makes
        // the code a lot cleaner.
        object_coordinates.resize(number_of_objects);
        temperature_submodules.resize(number_of_objects);
        temperatures.resize(number_of_objects);
        temperature_depths.resize(number_of_objects);
        composition_submodules.resize(number_of_objects);
        composition_depths.resize(number_of_objects);
        compositions.resize(number_of_objects);

        for (unsigned int object_number = 0; object_number < number_of_objects; object_number++)
          {
            wb_prm.enter_subsection("array[" + Utilities::to_string(object_number) + "]");
            {
              bool collect = false;
              wb_prm.enter_subsection("module");
              {
                collect = (boost::algorithm::trim_copy(boost::algorithm::to_lower_copy(wb_prm.get("name"))) == "continental plate");
              }
              wb_prm.leave_subsection();
              if (collect)
                {
                  object_coordinates[object_number] = Utilities::vector_double_to_point<2>(Utilities::vector_vector_string_to_double(wb_prm.get_double_array("coordinates")));

                  wb_prm.enter_subsection("module");
                  {
                    wb_prm.enter_subsection("temperature submodule");
                    {
                      temperature_submodules[object_number] = boost::algorithm::trim_copy(boost::algorithm::to_lower_copy(wb_prm.get("name")));
                      temperatures[object_number] = Utilities::string_to_double(wb_prm.get("temperature"));
                      temperature_depths[object_number] = Utilities::string_to_double(wb_prm.get("depth"));
                    }
                    wb_prm.leave_subsection();
                    wb_prm.enter_subsection("composition submodule");
                    {
                      composition_submodules[object_number] = boost::algorithm::trim_copy(boost::algorithm::to_lower_copy(wb_prm.get("name")));
                      compositions[object_number] = Utilities::string_to_double(wb_prm.get("composition"));
                      composition_depths[object_number] = Utilities::string_to_double(wb_prm.get("depth"));
                    }
                    wb_prm.leave_subsection();
                  }
                  wb_prm.leave_subsection();
                }
            }
            wb_prm.leave_subsection();
          }
      }
      wb_prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace WorldBuilder
  {
    ASPECT_REGISTER_WORLD_BUILDER(ContinentalPlate,
                                  "continental plate",
                                  "Implementation of a model in which the initial topography "
                                  "is zero. ")
  }
}
