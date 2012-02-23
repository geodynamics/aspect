//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/interface.h>
#include <aspect/simulator.h>

#include <typeinfo>


namespace aspect
{
  namespace Postprocess
  {
// ------------------------------ Interface -----------------------------

    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::declare_parameters (ParameterHandler &)
    {}



    template <int dim>
    void
    Interface<dim>::parse_parameters (ParameterHandler &)
    {}



    template <int dim>
    void
    Interface<dim>::save (std::map<std::string,std::string> &) const
    {}


    template <int dim>
    void
    Interface<dim>::load (const std::map<std::string,std::string> &)
    {}


// ------------------------------ SimulatorAccess -----------------------

    template <int dim>
    void
    SimulatorAccess<dim>::initialize (const Simulator<dim> &simulator_object)
    {
      simulator = SmartPointer<const Simulator<dim> > (&simulator_object, typeid(*this).name());
    }



    template <int dim>
    double SimulatorAccess<dim>::get_time () const
    {
      return simulator->time;
    }

    template <int dim>
    const ConditionalOStream &
    SimulatorAccess<dim>::get_pcout () const
    {
      return simulator->pcout;
    }

    template <int dim>
    double SimulatorAccess<dim>::get_timestep () const
    {
      return simulator->time_step;
    }



    template <int dim>
    unsigned int SimulatorAccess<dim>::get_timestep_number () const
    {
      return simulator->timestep_number;
    }



    template <int dim>
    const parallel::distributed::Triangulation<dim> &
    SimulatorAccess<dim>::get_triangulation () const
    {
      return simulator->triangulation;
    }



    template <int dim>
    double
    SimulatorAccess<dim>::get_volume () const
    {
      return simulator->global_volume;
    }



    template <int dim>
    const Mapping<dim> &
    SimulatorAccess<dim>::get_mapping () const
    {
      return simulator->mapping;
    }



    template <int dim>
    std::string
    SimulatorAccess<dim>::get_output_directory () const
    {
      return simulator->parameters.output_directory;
    }



    template <int dim>
    bool
    SimulatorAccess<dim>::convert_output_to_years () const
    {
      return simulator->parameters.convert_to_years;
    }



    template <int dim>
    void
    SimulatorAccess<dim>::get_refinement_criteria (Vector<float> & estimated_error_per_cell) const
    {
      simulator->compute_refinement_criterion(estimated_error_per_cell);
    }



    template <int dim>
    const TrilinosWrappers::MPI::BlockVector &
    SimulatorAccess<dim>::get_solution () const
    {
      return simulator->system_solution;
    }



    template <int dim>
    const TrilinosWrappers::MPI::BlockVector &
    SimulatorAccess<dim>::get_old_solution () const
    {
      return simulator->old_system_solution;
    }



    template <int dim>
    const DoFHandler<dim> &
    SimulatorAccess<dim>::get_dof_handler () const
    {
      return simulator->system_dof_handler;
    }



    template <int dim>
    void
    SimulatorAccess<dim>::get_depth_average_temperature(std::vector<double> & values) const
    {
      simulator->compute_depth_average_temperature(values);
    }


    template <int dim>
    void
    SimulatorAccess<dim>::get_Vs_anomaly(Vector<float> & values) const
    {
      simulator->compute_Vs_anomaly(values);
    }


    template <int dim>
    const MaterialModel::Interface<dim> &
    SimulatorAccess<dim>::get_material_model () const
    {
      return *simulator->material_model.get();
    }



    template <int dim>
    const BoundaryTemperature::Interface<dim> &
    SimulatorAccess<dim>::get_boundary_temperature () const
    {
      return *simulator->boundary_temperature.get();
    }



    template <int dim>
    const GeometryModel::Interface<dim> &
    SimulatorAccess<dim>::get_geometry_model () const
    {
      return *simulator->geometry_model.get();
    }

    template <int dim>
    const GravityModel::Interface<dim> &
    SimulatorAccess<dim>::get_gravity_model () const
    {
      return *simulator->gravity_model.get();
    }


    template <int dim>
    const AdiabaticConditions<dim> &
    SimulatorAccess<dim>::get_adiabatic_conditions () const
    {
      return *simulator->adiabatic_conditions.get();
    }

// ------------------------------ Manager -----------------------------

    template <int dim>
    void
    Manager<dim>::initialize (const Simulator<dim> &simulator)
    {
      std::list<std::pair<std::string,std::string> > output_list;
      for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        dynamic_cast<SimulatorAccess<dim>&>(**p).initialize (simulator);
    }



    template <int dim>
    std::list<std::pair<std::string,std::string> >
    Manager<dim>::execute (TableHandler &statistics)
    {
      // call the execute() functions of all postprocessor objects we have
      // here in turns
      std::list<std::pair<std::string,std::string> > output_list;
      for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        {
          // call the execute() function. if it produces any output
          // then add it to the list
          std::pair<std::string,std::string> output
            = (*p)->execute (statistics);

          if (output.first.size() + output.second.size() > 0)
            output_list.push_back (output);
        }

      return output_list;
    }


// -------------------------------- Deal with registering postprocessors and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      internal::Plugins::PluginList<Interface<deal_II_dimension> > registered_plugins;

      template <>
      std::list<internal::Plugins::PluginList<Interface<deal_II_dimension> >::PluginInfo> *
      internal::Plugins::PluginList<Interface<deal_II_dimension> >::plugins = 0;
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
          = registered_plugins.get_pattern_of_names (true);
        prm.declare_entry("List of postprocessors",
                          "all",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of postprocessor objects that should be run "
                          "at the end of each time step. Some of these postprocessors will "
                          "declare their own parameters which may, for example, include that "
                          "they will actually do something only every so many time steps or "
                          "years. Alternatively, the text 'all' indicates that all available "
                          "postprocessors should be run after each time step.\n\n"
                          "The following postprocessors are available:\n\n"
                          +
                          registered_plugins.get_description_string());
      }
      prm.leave_subsection();

      // now declare the parameters of each of the registered
      // postprocessors in turn
      registered_plugins.declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      Assert (registered_plugins.plugins != 0,
              ExcMessage ("No postprocessors registered!?"));

      // first find out which postprocessors are requested
      std::vector<std::string> postprocessor_names;
      prm.enter_subsection("Postprocess");
      {
        postprocessor_names
          = Utilities::split_string_list(prm.get("List of postprocessors"));
      }
      prm.leave_subsection();

      // see if 'all' was selected (or is part of the list). if so
      // simply replace the list with one that contains all names
      if (std::find (postprocessor_names.begin(),
                     postprocessor_names.end(),
                     "all") != postprocessor_names.end())
        {
          postprocessor_names.clear();
          for (std::list<internal::Plugins::PluginList<Interface<deal_II_dimension> >::PluginInfo>::const_iterator
               p = registered_plugins.plugins->begin();
               p != registered_plugins.plugins->end(); ++p)
            postprocessor_names.push_back (std_cxx1x::get<0>(*p));
        }

      // then go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int name=0; name<postprocessor_names.size(); ++name)
        postprocessors.push_back (std_cxx1x::shared_ptr<Interface<dim> >
                                  (registered_plugins.create_plugin (postprocessor_names[name],
                                                                     prm)));
    }


    template <int dim>
    void
    Manager<dim>::register_postprocessor (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          Interface<dim> * (*factory_function) ())
    {
      registered_plugins.register_plugin (name,
                                          description,
                                          declare_parameters_function,
                                          factory_function);
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class Interface<deal_II_dimension>;
    template class SimulatorAccess<deal_II_dimension>;
    template class Manager<deal_II_dimension>;

    Manager<deal_II_dimension> m;
  }
}
