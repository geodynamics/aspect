/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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
#include <aspect/initial_composition/world_builder.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/global.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <world_builder/world.h>

namespace WorldBuilder
{
  class World;
}


namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file through the
     * World builder.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class WorldBuilderDisequilibrium : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        /**
         * Return the initial composition as a function of position and index
         * of the compositional field. If c is equal to the index of the bound_fluid
         * compositional field, then this function will increase the value of the
         * bound_fluid field by a user defined percentage.
         */
        double initial_composition (const Point<dim> &position, const unsigned int c) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A vector that specifies for each compositional field index if the world builder
         * should be evaluated for this compositional field.
         */
        std::vector<bool> compositional_field_defined_by_world_builder;

        /**
         * A pointer to the WorldBuilder object. Keeping this pointer ensures
         * that the object doesn't go away while we still need it.
         */
        std::shared_ptr<const ::WorldBuilder::World> world_builder;

        /**
         * A value between 0 and 100, that defines the percent increase
         * in the bound_fluid from it's equilibrium value.
        */
        double percentage_increase;
    };
  }
}


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    WorldBuilderDisequilibrium<dim>::
    initialize()
    {
      CitationInfo::add("GWB");
      world_builder = this->get_world_builder_pointer();
    }



    template <int dim>
    double
    WorldBuilderDisequilibrium<dim>::
    initial_composition (const Point<dim> &position, const unsigned int c) const
    {
      const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");
      if (compositional_field_defined_by_world_builder[c] == true)
        {
          double comp_value = world_builder->composition(Utilities::convert_point_to_array(position),
                                                         -this->get_geometry_model().height_above_reference_surface(position),
                                                         c);
          if (c == bound_fluid_idx &&
              this->get_parameters().compositional_field_methods[bound_fluid_idx] == Parameters<dim>::AdvectionFieldMethod::static_field)
            comp_value *= (1 + percentage_increase/100);
          return comp_value;
        }

      return 0.0;
    }



    template <int dim>
    void
    WorldBuilderDisequilibrium<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("World builder");
        {
          prm.declare_entry ("List of relevant compositions", "",
                             Patterns::Anything(),
                             "A list of names of compositional fields for "
                             "which to determine the initial composition using "
                             "the World Builder. As World Builder evaluations can "
                             "be expensive, this parameter allows to only evaluate "
                             "the fields that are relevant. This plugin returns 0.0 "
                             "for all compositions that are not selected in the list. "
                             "By default the list is empty and the world builder is "
                             "evaluated for all compositional fields.");
          prm.declare_entry ("Disequilibrium percentage", "5",
                             Patterns::Double(0),
                             "A percentage value that defines the increase in the "
                             "bound_fluid compositional field from its equilibrium value. "
                             "This is only applied if the bound_fluid compositional field "
                             "is defined by the World Builder. The default value is 10%, "
                             "which means that the initial value of the bound_fluid field "
                             "is 10 percent larger than its equilibrium value.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    WorldBuilderDisequilibrium<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("World builder");
        {
          const std::vector<std::string> composition_names = Utilities::split_string_list(prm.get("List of relevant compositions"));
          percentage_increase = prm.get_double("Disequilibrium percentage");

          compositional_field_defined_by_world_builder.resize(this->n_compositional_fields(),false);

          if (composition_names.size() == 0)
            for (unsigned int i=0; i<this->n_compositional_fields(); ++i)
              compositional_field_defined_by_world_builder[i] = true;

          for (const auto &composition_name: composition_names)
            {
              AssertThrow(this->introspection().compositional_name_exists (composition_name),
                          ExcMessage("All fields in \"List of relevant compositions\" must match names of compositional "
                                     "fields as assigned in the \"Compositional fields/Names of fields\" parameter."));

              compositional_field_defined_by_world_builder[this->introspection().compositional_index_for_name(composition_name)] = true;
            }
        }

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(WorldBuilderDisequilibrium,
                                              "world builder disequilibrium",
                                              "Specify the initial composition through the World Builder. "
                                              "If the composition is not named bound_fluid, the composition "
                                              "is set to the value returned by the World Builder. If the "
                                              "composition is named bound_fluid, the value returned by the "
                                              "World Builder is increased by a user defined percentage. ")
  }
}
#endif
