//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/boundary_temperature/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>

    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &prm)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &prm)
    {}


// -------------------------------- Deal with registering models and automating
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
    register_boundary_temperature (const std::string &name,
                                   const std::string &description,
                                   void (*declare_parameters_function) (ParameterHandler &),
                                   Interface<dim> * (*factory_function) ())
    {
      // see if this is the first time we get into this
      // function and if so initialize the variable above
      if (registered_plugins.plugins == 0)
        registered_plugins.plugins = new std::list<internal::Plugins::PluginList<Interface<deal_II_dimension> >::PluginInfo>();

      // now add one record to the list
      registered_plugins.plugins->push_back (internal::Plugins::PluginList<Interface<deal_II_dimension> >::PluginInfo(name,
                                                                           description,
                                                                           declare_parameters_function,
                                                                           factory_function));
    }


    template <int dim>
    Interface<dim> *
    create_boundary_temperature (ParameterHandler &prm)
    {
      Assert (registered_plugins.plugins != 0, ExcInternalError());

      std::string model_name;
      prm.enter_subsection ("Boundary temperature model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      for (std::list<internal::Plugins::PluginList<Interface<deal_II_dimension> >::PluginInfo>::const_iterator p = registered_plugins.plugins->begin();
           p != registered_plugins.plugins->end(); ++p)
        if (std_cxx1x::get<0>(*p) == model_name)
          {
            Interface<dim> *i = std_cxx1x::get<3>(*p)();
            i->parse_parameters (prm);
            return i;
          }

      AssertThrow (false, ExcNotImplemented());
      return 0;
    }



    void
    declare_parameters (ParameterHandler &prm)
    {
      Assert (registered_plugins.plugins != 0, ExcInternalError());

      // first collect a list of all registered models
      std::string model_names;
      for (std::list<internal::Plugins::PluginList<Interface<deal_II_dimension> >::PluginInfo>::const_iterator p = registered_plugins.plugins->begin();
           p != registered_plugins.plugins->end(); ++p)
        {
          if (model_names.size() > 0)
            model_names += "|";
          model_names += std_cxx1x::get<0>(*p);
        }

      // then declare the actual entry in the parameter file
      prm.enter_subsection ("Boundary temperature model");
      {
        prm.declare_entry ("Model name", "",
                           Patterns::Selection (model_names),
                           "Select one of the available boundary temperature models.");
      }
      prm.leave_subsection ();

      for (std::list<internal::Plugins::PluginList<Interface<deal_II_dimension> >::PluginInfo>::const_iterator p = registered_plugins.plugins->begin();
           p != registered_plugins.plugins->end(); ++p)
        std_cxx1x::get<2>(*p)(prm);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    template class Interface<deal_II_dimension>;

    template
    void
    register_boundary_temperature<deal_II_dimension> (const std::string &,
                                                      const std::string &,
                                                      void ( *) (ParameterHandler &),
                                                      Interface<deal_II_dimension> * ( *) ());

    template
    Interface<deal_II_dimension> *
    create_boundary_temperature<deal_II_dimension> (ParameterHandler &prm);
  }
}
