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


#include <aspect/termination_criteria/interface.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <typeinfo>

namespace aspect
{
  namespace TerminationCriteria
  {
// ------------------------------ Interface -----------------------------

    template <int dim>
    Interface<dim>::~Interface ()
    {}

    template <int dim>
    void
    Interface<dim>::initialize ()
    {}

    template <int dim>
    void
    Interface<dim>::declare_parameters (ParameterHandler &)
    {}

    template <int dim>
    double Interface<dim>::check_for_last_time_step (const double time_step) const
    {
      return time_step;
    }

    template <int dim>
    void
    Interface<dim>::parse_parameters (ParameterHandler &)
    {}



// ------------------------------ Manager -----------------------------

    template <int dim>
    double Manager<dim>::check_for_last_time_step (const double time_step) const
    {
      double new_time_step = time_step;
      for (const auto &p : termination_objects)
        {
          double current_time_step = p->check_for_last_time_step (new_time_step);

          AssertThrow (current_time_step > 0,
                       ExcMessage("Time step must be greater than 0."));
          AssertThrow (current_time_step <= new_time_step,
                       ExcMessage("Current time step must be less than or equal to time step entered into function."));

          new_time_step = std::min(current_time_step, new_time_step);
        }
      return new_time_step;
    }

    template <int dim>
    bool
    Manager<dim>::execute () const
    {
      bool terminate_simulation = false;


      // call the execute() functions of all plugins we have
      // here in turns.
      std::list<std::string>::const_iterator  itn = termination_obj_names.begin();
      for (typename std::list<std::unique_ptr<Interface<dim> > >::const_iterator
           p = termination_objects.begin();
           p != termination_objects.end(); ++p,++itn)
        {
          try
            {
              const bool terminate = (*p)->execute ();

              // do the reduction: does any one of the processors
              // think that we should terminate? (do the reduction in
              // data type int since there is currently no function
              // Utilities::MPI::CollectiveOr or similar)
              const bool all_terminate = (Utilities::MPI::max ((terminate ? 1 : 0),
                                                               this->get_mpi_communicator())
                                          == 1);
              terminate_simulation |= all_terminate;

              // Let the user know which criterion caused the termination
              if (all_terminate == true)
                this->get_pcout() << "Termination requested by criterion: "
                                  << *itn
                                  << std::endl;
            }
          // plugins that throw exceptions usually do not result in
          // anything good because they result in an unwinding of the stack
          // and, if only one processor triggers an exception, the
          // destruction of objects often causes a deadlock. thus, if
          // an exception is generated, catch it, print an error message,
          // and abort the program
          catch (std::exception &exc)
            {
              std::cerr << std::endl << std::endl
                        << "----------------------------------------------------"
                        << std::endl;
              std::cerr << "Exception on MPI process <"
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << "> while running termination criterion plugin <"
                        << typeid(**p).name()
                        << ">: " << std::endl
                        << exc.what() << std::endl
                        << "Aborting!" << std::endl
                        << "----------------------------------------------------"
                        << std::endl;

              // terminate the program!
              MPI_Abort (MPI_COMM_WORLD, 1);
            }
          catch (...)
            {
              std::cerr << std::endl << std::endl
                        << "----------------------------------------------------"
                        << std::endl;
              std::cerr << "Exception on MPI process <"
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << "> while running termination criterion plugin <"
                        << typeid(**p).name()
                        << ">: " << std::endl;
              std::cerr << "Unknown exception!" << std::endl
                        << "Aborting!" << std::endl
                        << "----------------------------------------------------"
                        << std::endl;

              // terminate the program!
              MPI_Abort (MPI_COMM_WORLD, 1);
            }
        }

      return terminate_simulation;
    }


// -------------------------------- Deal with registering plugins and automating
// -------------------------------- their setup and selection at run time

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
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // first declare the postprocessors we know about to
      // choose from
      prm.enter_subsection("Termination criteria");
      {
        // construct a string for Patterns::MultipleSelection that
        // contains the names of all registered termination criteria
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();
        prm.declare_entry("Termination criteria",
                          "end time",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of termination criteria that will "
                          "determine when the simulation should end. "
                          "Whether explicitly stated or not, the ``end time'' "
                          "termination criterion will always be used."
                          "The following termination criteria are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection();

      // now declare the parameters of each of the registered
      // plugins in turn
      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      Assert (std::get<dim>(registered_plugins).plugins != nullptr,
              ExcMessage ("No termination criteria plugins registered!?"));

      // first find out which plugins are requested
      std::vector<std::string> plugin_names;
      prm.enter_subsection("Termination criteria");
      {
        plugin_names = Utilities::split_string_list(prm.get("Termination criteria"));
        AssertThrow(Utilities::has_unique_entries(plugin_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Termination criteria/Termination criteria' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        // as described, the end time plugin is always active
        if (std::find (plugin_names.begin(), plugin_names.end(), "end time")
            == plugin_names.end())
          plugin_names.emplace_back("end time");
      }
      prm.leave_subsection();

      // go through the list, create objects, initialize them, and let them parse
      // their own parameters
      for (unsigned int name=0; name<plugin_names.size(); ++name)
        {
          termination_objects.push_back (std::unique_ptr<Interface<dim> >
                                         (std::get<dim>(registered_plugins)
                                          .create_plugin (plugin_names[name],
                                                          "Termination criteria::Termination criteria")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*termination_objects.back()))
            sim->initialize_simulator (this->get_simulator());
          termination_objects.back()->parse_parameters (prm);
          termination_objects.back()->initialize ();

          termination_obj_names.push_back(plugin_names[name]);
        }
    }


    template <int dim>
    void
    Manager<dim>::register_termination_criterion (const std::string &name,
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
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Termination criteria interface",
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
      std::list<internal::Plugins::PluginList<TerminationCriteria::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<TerminationCriteria::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<TerminationCriteria::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<TerminationCriteria::Interface<3> >::plugins = nullptr;
    }
  }

  namespace TerminationCriteria
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
