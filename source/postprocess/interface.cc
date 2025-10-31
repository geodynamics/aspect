/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/interface.h>
#include <aspect/utilities.h>

#include <typeinfo>


namespace aspect
{
  namespace Postprocess
  {
// ------------------------------ Interface -----------------------------

    template <int dim>
    std::list<std::string>
    Interface<dim>::required_other_postprocessors() const
    {
      return {};
    }



// ------------------------------ Manager -----------------------------


    template <int dim>
    std::list<std::pair<std::string,std::string>>
    Manager<dim>::execute (TableHandler &statistics)
    {
      // call the execute() functions of all postprocessor objects we have
      // here in turns
      std::list<std::pair<std::string,std::string>> output_list;
      for (const auto &p : this->plugin_objects)
        {
          try
            {
              // first call the update() function.
              p->update();

              // call the execute() function. if it produces any output
              // then add it to the list
              std::pair<std::string,std::string> output
                = p->execute (statistics);

              if (output.first.size() + output.second.size() > 0)
                output_list.push_back (output);
            }
          // postprocessors that throw exceptions usually do not result in
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
                        << "> while running postprocessor <"
                        << typeid(*p).name()
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
                        << "> while running postprocessor <"
                        << typeid(*p).name()
                        << ">: " << std::endl;
              std::cerr << "Unknown exception!" << std::endl
                        << "Aborting!" << std::endl
                        << "----------------------------------------------------"
                        << std::endl;

              // terminate the program!
              MPI_Abort (MPI_COMM_WORLD, 1);
            }
        }

      return  output_list;
    }


// -------------------------------- Deal with registering postprocessors and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::PluginList<Interface<2>>,
      aspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // first declare the postprocessors we know about to
      // choose from
      prm.enter_subsection("Postprocess");
      {
        // construct a string for Patterns::MultipleSelection that
        // contains the names of all registered postprocessors
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();
        prm.declare_entry("List of postprocessors",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of postprocessor objects that should be run "
                          "at the end of each time step. Some of these postprocessors will "
                          "declare their own parameters which may, for example, include that "
                          "they will actually do something only every so many time steps or "
                          "years. Alternatively, the text `all' indicates that all available "
                          "postprocessors should be run after each time step.\n\n"
                          "The following postprocessors are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection();

      // now declare the parameters of each of the registered
      // postprocessors in turn
      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      Assert (std::get<dim>(registered_plugins).plugins != nullptr,
              ExcMessage ("No postprocessors registered!?"));

      // first find out which postprocessors are requested
      prm.enter_subsection("Postprocess");
      {
        this->plugin_names
          = Utilities::split_string_list(prm.get("List of postprocessors"));
        AssertThrow(Utilities::has_unique_entries(this->plugin_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Postprocess/List of postprocessors' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));
      }
      prm.leave_subsection();

      // see if 'all' was selected (or is part of the list). if so
      // simply replace the list with one that contains all names
      if (std::find (this->plugin_names.begin(),
                     this->plugin_names.end(),
                     "all") != this->plugin_names.end())
        {
          this->plugin_names.clear();
          for (const auto &p : *std::get<dim>(registered_plugins).plugins)
            this->plugin_names.push_back (std::get<0>(p));
        }

      // see if the user specified "global statistics" somewhere; if so, remove
      // it from the list because we will *always* want to have it and so
      // whether or not it has been explicitly provided by the user makes no
      // difference.
      std::vector<std::string>::iterator new_end
        = std::remove (this->plugin_names.begin(),
                       this->plugin_names.end(),
                       "global statistics");
      if (new_end != this->plugin_names.end())
        this->plugin_names.erase (new_end, this->plugin_names.end());

      // in any case, put the global statistics postprocessor at the front:
      this->plugin_names.insert(this->plugin_names.begin(), "global statistics");

      // then go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int name=0; name<this->plugin_names.size(); ++name)
        {
          this->plugin_objects.emplace_back (std::get<dim>(registered_plugins)
                                             .create_plugin (this->plugin_names[name],
                                                             "Postprocessor plugins"));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*this->plugin_objects.back()))
            sim->initialize_simulator (this->get_simulator());

          this->plugin_objects.back()->parse_parameters (prm);
          this->plugin_objects.back()->initialize ();

          // now see if the newly created postprocessor relies on others. if so,
          // go through the list of the ones we already have and if the required
          // ones are new, add them to the end of the list we work through
          const std::list<std::string> additional_postprocessors
            = this->plugin_objects.back()->required_other_postprocessors ();

          for (const auto &p : additional_postprocessors)
            {
              AssertThrow (Patterns::Selection(std::get<dim>(registered_plugins).get_pattern_of_names ())
                           .match (p) == true,
                           ExcMessage ("Postprocessor <" + this->plugin_names[name] +
                                       "> states that it depends on another postprocessor, <"
                                       + p +
                                       ">, but the latter is not a valid name."));

              bool already_present = false;
              for (const auto &postprocessor_name : this->plugin_names)
                if (postprocessor_name == p)
                  {
                    already_present = true;
                    break;
                  }

              if (already_present == false)
                this->plugin_names.push_back (p);
            }
        }
      Assert (this->plugin_names.size() == this->plugin_objects.size(),
              ExcInternalError());

      // we now have matching lists 'this->plugin_objects' and 'this->plugin_names'
      // that define which postprocessors we have. we just need to bring them
      // into an order so that dependencies run first, and dependents after
      // them. the algorithm for creating a sorted list for this works as follows:
      //
      // while there are postprocessors not yet added to the sorted list:
      // - go through the list
      // - if we encounter a postprocessor that has not yet been added yet:
      //   . if all of its dependencies are in the list, add it to the end
      //   . if at least one of its dependencies are not yet in the list,
      //     skip it
      //
      // if we go through a loop where we do not add a postprocessor to the list
      // but there are still ones that haven't been added to the list, then we
      // have found a cycle in the dependencies and that is clearly a problem
      std::vector<bool> already_assigned (this->plugin_objects.size(), false);
      std::vector<std::string> sorted_names;
      std::list<std::unique_ptr<Interface<dim>>> sorted_postprocessors;
      while (sorted_names.size() < this->plugin_objects.size())
        {
          bool at_least_one_element_added = false;

          {
            typename std::list<std::unique_ptr<Interface<dim>>>::iterator
            pp = this->plugin_objects.begin();
            for (unsigned int i=0; i<this->plugin_names.size(); ++i, ++pp)
              if (already_assigned[i] == false)
                {
                  // for this postprocessor, check if all of its dependencies
                  // are already in the list
                  const std::list<std::string> deps = (*pp)->required_other_postprocessors();
                  bool unmet_dependencies = false;
                  for (const auto &p : deps)
                    if (std::find (sorted_names.begin(),
                                   sorted_names.end(),
                                   p) == sorted_names.end())
                      {
                        unmet_dependencies = true;
                        break;
                      }

                  // if we have unmet dependencies, there is nothing we can do
                  // right now for this postprocessor (but we will come back for it)
                  //
                  // if there are none, add this postprocessor. using move semantics
                  // removes the postprocessor from the 'this->plugin_objects' array,
                  // which is ok since we will swap the two arrays at the end of the
                  // function.
                  if (unmet_dependencies == false)
                    {
                      sorted_names.push_back (this->plugin_names[i]);
                      sorted_postprocessors.emplace_back (std::move(*pp));
                      already_assigned[i] = true;
                      at_least_one_element_added = true;
                    }
                }
          }

          // check that we have added at least one element; if not, there is a cycle
          if (at_least_one_element_added == false)
            {
              std::ostringstream out;
              out << "While sorting postprocessors by their dependencies, "
                  "ASPECT encountered a cycle in dependencies. The following "
                  "postprocessors are involved:\n";
              typename std::list<std::unique_ptr<Interface<dim>>>::const_iterator
              pp = this->plugin_objects.begin();
              for (unsigned int i=0; i<this->plugin_names.size(); ++i, ++pp)
                if (already_assigned[i] == false)
                  {
                    out << "  " << this->plugin_names[i] << " -> ";
                    const std::list<std::string> deps = (*pp)->required_other_postprocessors();
                    for (const auto &p : deps)
                      out << "'" << p << "' ";
                    out << std::endl;
                  }
              AssertThrow (false, ExcMessage(out.str()));
            }
        }
      Assert (this->plugin_names.size() == sorted_names.size(),
              ExcInternalError());
      Assert (sorted_postprocessors.size() == sorted_names.size(),
              ExcInternalError());
      Assert (std::find (already_assigned.begin(), already_assigned.end(), false) == already_assigned.end(),
              ExcInternalError());

      // finally swap the unsorted lists with the sorted lists and only
      // keep the latter
      this->plugin_objects.swap (sorted_postprocessors);
      this->plugin_names.swap (sorted_names);
    }


    template <int dim>
    void
    Manager<dim>::register_postprocessor (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          std::unique_ptr<Interface<dim>> (*factory_function) ())
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
      std::get<dim>(registered_plugins).write_plugin_graph ("Postprocessor interface",
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
      std::list<internal::Plugins::PluginList<Postprocess::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<Postprocess::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::Interface<3>>::plugins = nullptr;
    }
  }

  namespace Postprocess
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
