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

#include <aspect/global.h>

#ifdef ASPECT_WITH_WORLD_BUILDER
#  include <aspect/geometry_model/interface.h>
#  include <aspect/prescribed_solution/world_builder.h>

#  include <world_builder/world.h>

#  include <algorithm>

namespace aspect
{
  namespace PrescribedSolution
  {
    namespace
    {
      constexpr unsigned int temperature_indicator = 0;
      constexpr unsigned int velocity_indicator = 1;
      constexpr unsigned int composition_indicator = 2;



      template <typename WorldBuilderType>
      auto
      supports_indicators (const WorldBuilderType &world_builder, int)
      -> decltype(world_builder.indicator(std::array<double,2>(), 0.0, 0), bool())
      {
        return true;
      }



      template <typename WorldBuilderType>
      bool
      supports_indicators (const WorldBuilderType &, long)
      {
        return false;
      }



      template <typename WorldBuilderType, typename PositionType>
      auto
      get_indicator (const WorldBuilderType &world_builder,
                     const PositionType &position,
                     const double depth,
                     const unsigned int indicator,
                     int)
      -> decltype(world_builder.indicator(position, depth, indicator))
      {
        return world_builder.indicator(position, depth, indicator);
      }



      template <typename WorldBuilderType, typename PositionType>
      double
      get_indicator (const WorldBuilderType &,
                     const PositionType &,
                     const double,
                     const unsigned int,
                     long)
      {
        AssertThrow(false,
                    ExcMessage("The prescribed solution model 'world builder' requires "
                               "a WorldBuilder version that supports indicator fields."));
        return 0.0;
      }
    }



    template <int dim>
    void
    WorldBuilder<dim>::initialize ()
    {
      CitationInfo::add("GWB");
      world_builder = this->get_world_builder_pointer();

      AssertThrow(supports_indicators(*world_builder, 0),
                  ExcMessage("The prescribed solution model 'world builder' requires "
                             "a WorldBuilder version that supports indicator fields."));
    }



    template <int dim>
    void
    WorldBuilder<dim>::declare_parameters (ParameterHandler &)
    {}



    template <int dim>
    void
    WorldBuilder<dim>::constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &,
                                           const std::vector<Point<dim>> &positions,
                                           const std::vector<unsigned int> &component_indices,
                                           std::vector<bool> &should_be_constrained,
                                           std::vector<double> &solution)
    {
      const unsigned int temperature_component = this->introspection().component_indices.temperature;
      const unsigned int first_velocity_component = this->introspection().component_indices.velocities[0];
      const unsigned int last_velocity_component = this->introspection().component_indices.velocities[dim-1];
      const std::vector<unsigned int> &composition_components =
        this->introspection().component_indices.compositional_fields;

      for (unsigned int q=0; q<positions.size(); ++q)
        {
          const unsigned int component = component_indices[q];
          const std::array<double,dim> position = Utilities::convert_point_to_array(positions[q]);
          const double depth = -this->get_geometry_model().height_above_reference_surface(positions[q]);

          if (component == temperature_component)
            {
              if (get_indicator(*world_builder, position, depth, temperature_indicator, 0) > 0.5)
                {
                  solution[q] = world_builder->temperature(position, depth);
                  should_be_constrained[q] = true;
                }
            }
          else if (component >= first_velocity_component && component <= last_velocity_component)
            {
              if (get_indicator(*world_builder, position, depth, velocity_indicator, 0) > 0.5)
                {
                  const unsigned int direction = component - first_velocity_component;
                  const std::vector<double> velocity =
                  world_builder->properties(position, depth, {{{5,0,0}}});

                  solution[q] = velocity[direction] / year_in_seconds;
                  should_be_constrained[q] = true;
                }
            }
          else
            {
              const auto composition_component =
                std::find(composition_components.begin(), composition_components.end(), component);

              if (composition_component != composition_components.end() &&
                  get_indicator(*world_builder, position, depth, composition_indicator, 0) > 0.5)
                {
                  const unsigned int composition_index =
                    std::distance(composition_components.begin(), composition_component);

                  solution[q] = world_builder->composition(position, depth, composition_index);
                  should_be_constrained[q] = true;
                }
            }
        }
    }
  }
}



namespace aspect
{
  namespace PrescribedSolution
  {
    ASPECT_REGISTER_PRESCRIBED_SOLUTION(WorldBuilder,
                                        "world builder",
                                        "Prescribe velocity, temperature, and composition using values from "
                                        "the WorldBuilder. WorldBuilder indicators with indices 0, 1, and 2 "
                                        "select where temperature, velocity, and composition are constrained, "
                                        "respectively.")
  }
}

#endif
