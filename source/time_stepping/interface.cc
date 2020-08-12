/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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


#include <aspect/time_stepping/interface.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace TimeStepping
  {
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
    Interface<dim>::
    initialize ()
    {
    }



    template <int dim>
    void
    Interface<dim>::
    update ()
    {
    }



    template <int dim>
    void
    Interface<dim>::declare_parameters (ParameterHandler &/*prm*/)
    {
    }



    template <int dim>
    void
    Interface<dim>::
    parse_parameters (ParameterHandler &/*prm*/)
    {
    }



    template <int dim>
    double
    Manager<dim>::
    compute_time_step_size() const
    {
      double new_time_step = std::numeric_limits<double>::max();
      Reaction reaction = Reaction::advance;

      for (const auto &plugin : active_plugins)
        {
          std::pair<Reaction, double> answer
            = plugin->execute();

          new_time_step = std::min(new_time_step, answer.second);
          reaction = static_cast<Reaction>(std::min(reaction, answer.first));
        }

      // For now, we only support the default:
      AssertThrow(reaction == Reaction::advance,
                  ExcNotImplemented());

      // Make sure we do not exceed the maximum time step length. This can happen
      // if velocities get too small or even zero in models that are stably stratified
      // or use prescribed velocities.
      new_time_step = std::min(new_time_step, this->get_parameters().maximum_time_step);

      // Make sure that the time step doesn't increase too fast
      if (this->get_timestep() != 0)
        new_time_step = std::min(new_time_step, this->get_timestep() + this->get_timestep() * this->get_parameters().maximum_relative_increase_time_step);

      // Make sure we do not exceed the maximum length for the first time step
      if (this->get_timestep_number() == 0)
        new_time_step = std::min(new_time_step, this->get_parameters().maximum_first_time_step);

      // Make sure we reduce the time step length appropriately if we terminate after this step
      new_time_step = termination_manager.check_for_last_time_step(new_time_step);

      AssertThrow (new_time_step > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length was: " + std::to_string(new_time_step) + ". "
                              "Please check the time stepping plugins and for non-positive material properties."));
      return new_time_step;
    }



    template <int dim>
    bool
    Manager<dim>::
    need_checkpoint_on_terminate() const
    {
      return do_checkpoint_on_terminate;
    }



    template <int dim>
    bool
    Manager<dim>::should_simulation_terminate_now() const
    {
      return termination_manager.execute();
    }



    template <int dim>
    void
    Manager<dim>::initialize_simulator (const Simulator<dim> &simulator_object)
    {
      SimulatorAccess<dim>::initialize_simulator(simulator_object);
      termination_manager.initialize_simulator(simulator_object);
    }



    template <int dim>
    void
    Manager<dim>::register_time_stepping_model(const std::string &name,
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
    Manager<dim>:: declare_parameters (ParameterHandler &prm)
    {
      TerminationCriteria::Manager<dim>::declare_parameters(prm);
      prm.enter_subsection("Termination criteria");
      {
        prm.declare_entry("Checkpoint on termination", "false",
                          Patterns::Bool (),
                          "Whether to checkpoint the simulation right before termination.");
      }
      prm.leave_subsection();


      prm.enter_subsection("Time stepping");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of time stepping plugins that "
                          "will be used to calculate the time step size. The minimum of the "
                          " result of each plugin will be used.\n\n"
                          "The following plugins are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());

      }
      prm.leave_subsection();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      termination_manager.parse_parameters(prm);
      prm.enter_subsection("Termination criteria");
      {
        do_checkpoint_on_terminate = prm.get_bool("Checkpoint on termination");
      }
      prm.leave_subsection();

      {
        prm.enter_subsection("Time stepping");
        std::vector<std::string>
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Time stepping/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        if (model_names.size()==0)
          {
            // handle the default case, where the user has not chosen any time stepping scheme explicitly:

            model_names.emplace_back("convection time step");
            if (this->get_parameters().use_conduction_timestep)
              model_names.emplace_back("conduction time step");
          }
        else
          {
            AssertThrow(this->get_parameters().use_conduction_timestep == false,
                        ExcMessage("When you are using Time stepping:: List of model names, do can not "
                                   "set \"Use conduction timestep\" to true. Use the \"conduction time step\" "
                                   "instead"));
          }

        prm.leave_subsection();

        for (unsigned int name=0; name<model_names.size(); ++name)
          {
            active_plugins.push_back (std::unique_ptr<Interface<dim> >
                                      (std::get<dim>(registered_plugins)
                                       .create_plugin (model_names[name],
                                                       "Time stepping::Model names")));

            if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*active_plugins.back()))
              sim->initialize_simulator (this->get_simulator());

            active_plugins.back()->parse_parameters (prm);
            active_plugins.back()->initialize ();
          }
      }
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Time stepping interface",
                                                            out);
    }


    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }

  }
}




// explicit instantiation of the functions we implement in this file
namespace aspect
{

  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<TimeStepping::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<TimeStepping::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<TimeStepping::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<TimeStepping::Interface<3> >::plugins = nullptr;
    }
  }

  namespace TimeStepping
  {
#define INSTANTIATE(dim) \
  \
  template class Interface<dim>; \
  template class Manager<dim>; \
  template \
  std::string \
  get_valid_model_names_pattern<dim> ();

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
