//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/gravity_model/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace GravityModel
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


// -------------------------------- Deal with registering gravity models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      internal::Plugins::PluginList<Interface<deal_II_dimension> > registered_plugins;
    }



    template <int dim>
    void
    register_gravity_model (const std::string &name,
                            const std::string &description,
                            void (*declare_parameters_function) (ParameterHandler &),
                            Interface<dim> *(*factory_function) ())
    {
      registered_plugins.register_plugin (name,
                                          description,
                                          declare_parameters_function,
                                          factory_function);
    }


    template <int dim>
    Interface<dim> *
    create_gravity_model (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Gravity model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      return registered_plugins.create_plugin (model_name, prm);
    }



    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Gravity model");
      {
        const std::string pattern_of_names
          = registered_plugins.get_pattern_of_names ();
        prm.declare_entry ("Model name", "",
                           Patterns::Selection (pattern_of_names),
                           "Select one of the following models:\n\n"
                           +
                           registered_plugins.get_description_string());
      }
      prm.leave_subsection ();

      registered_plugins.declare_parameters (prm);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace internal {
    namespace Plugins {
      template <>
      std::list<internal::Plugins::PluginList<GravityModel::Interface<deal_II_dimension> >::PluginInfo> *
      internal::Plugins::PluginList<GravityModel::Interface<deal_II_dimension> >::plugins = 0;
    }
  }

  namespace GravityModel
  {
    template class Interface<deal_II_dimension>;

    template
    void
    register_gravity_model<deal_II_dimension> (const std::string &,
                                               const std::string &,
                                               void ( *) (ParameterHandler &),
                                               Interface<deal_II_dimension> *( *) ());

    template
    Interface<deal_II_dimension> *
    create_gravity_model<deal_II_dimension> (ParameterHandler &prm);
  }
}
