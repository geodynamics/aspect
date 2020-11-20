/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
#include <aspect/boundary_composition/interface.h>

#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>
#include <tuple>

#include <list>


namespace aspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}

    template <int dim>
    void
    Interface<dim>::update ()
    {}

    template <int dim>
    void
    Interface<dim>::initialize ()
    {}

    template <int dim>
    double
    Interface<dim>::boundary_composition (const types::boundary_id /*boundary_indicator*/,
                                          const Point<dim>        &/*position*/,
                                          const unsigned int       /*compositional_field*/) const
    {
      AssertThrow(false,
                  ExcMessage("The boundary composition plugin has to implement a function called `composition' "
                             "with four arguments or a function `boundary_composition' with three arguments. "
                             "The function with four arguments is deprecated and will "
                             "be removed in a later version of ASPECT."));
      return numbers::signaling_nan<double>();
    }

    template <int dim>
    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &)
    {}




    // ------------------------------ Manager -----------------------------
    // -------------------------------- Deal with registering boundary_composition models and automating
    // -------------------------------- their setup and selection at run time

    template <int dim>
    Manager<dim>::~Manager()
    {}



    template <int dim>
    void
    Manager<dim>::update ()
    {
      for (unsigned int i=0; i<boundary_composition_objects.size(); ++i)
        {
          boundary_composition_objects[i]->update();
        }
    }



    namespace
    {
      std::tuple
      <void *,
      void *,
      aspect::internal::Plugins::PluginList<Interface<2> >,
      aspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_boundary_composition (const std::string &name,
                                                 const std::string &description,
                                                 void (*declare_parameters_function) (ParameterHandler &),
                                                 Interface<dim> *(*factory_function) ())
    {
      std::get<dim>(registered_plugins).register_plugin (name,
                                                         description,
                                                         declare_parameters_function,
                                                         factory_function);
    }


    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      // find out which plugins are requested and the various other
      // parameters we declare here
      prm.enter_subsection ("Boundary composition model");
      {
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Boundary composition model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        const std::string model_name = prm.get ("Model name");

        AssertThrow (model_name == "unspecified" || model_names.size() == 0,
                     ExcMessage ("The parameter 'Model name' is only used for reasons"
                                 "of backwards compatibility and can not be used together with "
                                 "the new functionality 'List of model names'. Please add your "
                                 "boundary composition model to the list instead."));

        if (!(model_name == "unspecified"))
          model_names.push_back(model_name);

        // create operator list
        std::vector<std::string> model_operator_names =
          Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list(prm.get("List of model operators")),
                                                  model_names.size(),
                                                  "List of model operators");
        model_operators = Utilities::create_model_operator_list(model_operator_names);

        try
          {
            const std::vector<types::boundary_id> x_fixed_composition_boundary_indicators
              = this->get_geometry_model().translate_symbolic_boundary_names_to_ids (Utilities::split_string_list
                                                                                     (prm.get ("Fixed composition boundary indicators")));
            fixed_composition_boundary_indicators
              = std::set<types::boundary_id> (x_fixed_composition_boundary_indicators.begin(),
                                              x_fixed_composition_boundary_indicators.end());

            // If model names have been set, but no boundaries on which to use them,
            // ignore the set values, do not create objects that are never used.
            if (fixed_composition_boundary_indicators.size() == 0)
              {
                model_names.clear();
                model_operators.clear();
              }
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Fixed composition "
                                            "boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows: "
                                            + error));
          }

        if (prm.get ("Allow fixed composition on outflow boundaries") == "true")
          allow_fixed_composition_on_outflow_boundaries = true;
        else if (prm.get ("Allow fixed composition on outflow boundaries") == "false")
          allow_fixed_composition_on_outflow_boundaries = false;
        else if (prm.get ("Allow fixed composition on outflow boundaries") == "false for models without melt")
          allow_fixed_composition_on_outflow_boundaries = this->get_parameters().include_melt_transport;
        else
          AssertThrow(false, ExcMessage("'Allow fixed composition on outflow boundaries' "
                                        "must be set to 'true' or 'false', or to its default value."));
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : model_names)
        {
          // create boundary composition objects
          boundary_composition_objects.push_back (std::unique_ptr<Interface<dim> >
                                                  (std::get<dim>(registered_plugins)
                                                   .create_plugin (model_name,
                                                                   "Boundary composition::Model names")));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_composition_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          boundary_composition_objects.back()->parse_parameters (prm);
          boundary_composition_objects.back()->initialize ();
        }
    }



    template <int dim>
    double
    Manager<dim>::boundary_composition (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position,
                                        const unsigned int compositional_field) const
    {
      double composition = 0.0;

      for (unsigned int i=0; i<boundary_composition_objects.size(); ++i)
        composition = model_operators[i](composition,
                                         boundary_composition_objects[i]->boundary_composition(boundary_indicator,
                                             position,
                                             compositional_field));

      return composition;
    }



    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_boundary_composition_names () const
    {
      return model_names;
    }


    template <int dim>
    const std::vector<std::unique_ptr<Interface<dim> > > &
    Manager<dim>::get_active_boundary_composition_conditions () const
    {
      return boundary_composition_objects;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_fixed_composition_boundary_indicators() const
    {
      return fixed_composition_boundary_indicators;
    }



    template <int dim>
    bool
    Manager<dim>::allows_fixed_composition_on_outflow_boundaries() const
    {
      return allow_fixed_composition_on_outflow_boundaries;
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Boundary composition model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of boundary composition models that "
                          "will be used to initialize the composition. "
                          "These plugins are loaded in the order given, and modify the "
                          "existing composition field via the operators listed "
                          "in 'List of model operators'.\n\n"
                          "The following boundary composition models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection(Utilities::get_model_operator_options()),
                          "A comma-separated list of operators that "
                          "will be used to append the listed composition models onto "
                          "the previous models. If only one operator is given, "
                          "the same operator is applied to all models.");

        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string()
                           + "\n\n" +
                           "\\textbf{Warning}: This parameter provides an old and "
                           "deprecated way of specifying "
                           "boundary composition models and shouldn't be used. "
                           "Please use 'List of model names' instead.");

        prm.declare_entry ("Fixed composition boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the composition is fixed and described by the "
                           "boundary composition object selected in its own section "
                           "of this input file. All boundary indicators used by the geometry "
                           "but not explicitly listed here will end up with no-flux "
                           "(insulating) boundary conditions."
                           "\n\n"
                           "The names of the boundaries listed here can either be "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model."
                           "\n\n"
                           "This parameter only describes which boundaries have a fixed "
                           "composition, but not what composition should hold on these "
                           "boundaries. The latter piece of information needs to be "
                           "implemented in a plugin in the BoundaryComposition "
                           "group, unless an existing implementation in this group "
                           "already provides what you want.");
        prm.declare_entry ("Allow fixed composition on outflow boundaries", "false for models without melt",
                           Patterns::Selection("true|false|false for models without melt"),
                           "When the composition is fixed on a given boundary as determined "
                           "by the list of 'Fixed composition boundary indicators', there "
                           "might be parts of the boundary where material flows out and "
                           "one may want to prescribe the composition only on those parts of "
                           "the boundary where there is inflow. This parameter determines "
                           "if compositions are only prescribed at these inflow parts of the "
                           "boundary (if false) or everywhere on a given boundary, independent "
                           "of the flow direction (if true). By default, this parameter is set "
                           "to false, except in models with melt transport (see below). "
                           "Note that in this context, `fixed' refers to the fact that these "
                           "are the boundary indicators where Dirichlet boundary conditions are "
                           "applied, and does not imply that the boundary composition is "
                           "time-independent. "
                           "\n\n"
                           "Mathematically speaking, the compositional fields satisfy an "
                           "advection equation that has no diffusion. For this equation, one "
                           "can only impose Dirichlet boundary conditions (i.e., prescribe a "
                           "fixed compositional field value at the boundary) at those boundaries "
                           "where material flows in. This would correspond to the ``false'' "
                           "setting of this parameter, which is correspondingly the default. "
                           "On the other hand, on a finite dimensional discretization such as "
                           "the one one obtains from the finite element method, it is possible "
                           "to also prescribe values on outflow boundaries, even though this may "
                           "make no physical sense. This would then correspond to the ``true'' "
                           "setting of this parameter."
                           "\n\n"
                           "A warning for models with melt transport: In models with fluid flow, "
                           "some compositional fields (in particular the porosity) might be "
                           "transported with the fluid velocity, and would need to set the "
                           "constraints based on the fluid velocity. However, this is currently "
                           "not possible, because we reuse the same matrix for all compositional "
                           "fields, and therefore can not use different constraints for different "
                           "fields. Consequently, we set this parameter to true by default in "
                           "models where melt transport is enabled. Be aware that if you change "
                           "this default setting, you will not use the melt velocity, but the solid "
                           "velocity to determine on which parts of the boundaries there is outflow.");
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Boundary composition interface",
                                                            out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<BoundaryComposition::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryComposition::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryComposition::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryComposition::Interface<3> >::plugins = nullptr;
    }
  }

  namespace BoundaryComposition
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
