//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/initial_conditions/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::initialize (const GeometryModel::Interface<dim>       &geometry_model_,
                                const BoundaryTemperature::Interface<dim> &boundary_temperature_,
                                const AdiabaticConditions<dim>            &adiabatic_conditions_)
    {
      geometry_model       = &geometry_model_;
      boundary_temperature = &boundary_temperature_;
      adiabatic_conditions = &adiabatic_conditions_;
    }


    template <int dim>
    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &prm)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &prm)
    {}


// -------------------------------- Deal with registering initial_conditions models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      typedef
      std_cxx1x::tuple<std::string,
                void ( *) (ParameterHandler &),
                Interface<deal_II_dimension> * ( *) ()>
                InitialConditionsModelInfo;

      // A pointer to a list of all postprocessors. the three elements of the tuple
      // correspond to the arguments given to the register_initial_conditions_model
      // function.
      //
      // The object is a pointer rather for the same reason as discussed in
      // postprocess_base.cc for the corresponding variable there
      std::list<InitialConditionsModelInfo> *registered_initial_conditions_models = 0;
    }



    template <int dim>
    void
    register_initial_conditions_model (const std::string &name,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> * (*factory_function) ())
    {
      // see if this is the first time we get into this
      // function and if so initialize the variable above
      if (registered_initial_conditions_models == 0)
        registered_initial_conditions_models = new std::list<InitialConditionsModelInfo>();

      // now add one record to the list
      registered_initial_conditions_models->push_back (InitialConditionsModelInfo(name,
                                                                                  declare_parameters_function,
                                                                                  factory_function));
    }


    template <int dim>
    Interface<dim> *
    create_initial_conditions (ParameterHandler &prm,
                               const GeometryModel::Interface<dim> &geometry_model,
                               const BoundaryTemperature::Interface<dim> &boundary_temperature,
                               const AdiabaticConditions<dim>      &adiabatic_conditions)
    {
      Assert (registered_initial_conditions_models != 0, ExcInternalError());

      std::string model_name;
      prm.enter_subsection ("Initial conditions");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      for (std::list<InitialConditionsModelInfo>::const_iterator p = registered_initial_conditions_models->begin();
           p != registered_initial_conditions_models->end(); ++p)
        if (std_cxx1x::get<0>(*p) == model_name)
          {
            Interface<dim> *i = std_cxx1x::get<2>(*p)();
            i->parse_parameters (prm);
            i->initialize (geometry_model,
                           boundary_temperature,
                           adiabatic_conditions);
            return i;
          }

      AssertThrow (false, ExcNotImplemented());
      return 0;
    }



    void
    declare_parameters (ParameterHandler &prm)
    {
      Assert (registered_initial_conditions_models != 0, ExcInternalError());

      // first collect a list of all registered models
      std::string model_names;
      for (std::list<InitialConditionsModelInfo>::const_iterator p = registered_initial_conditions_models->begin();
           p != registered_initial_conditions_models->end(); ++p)
        {
          if (model_names.size() > 0)
            model_names += "|";
          model_names += std_cxx1x::get<0>(*p);
        }

      // then declare the actual entry in the parameter file
      prm.enter_subsection ("Initial conditions");
      {
        prm.declare_entry ("Model name", "",
                           Patterns::Selection (model_names),
                           "Select one of the available initial conditions.");
      }
      prm.leave_subsection ();

      for (std::list<InitialConditionsModelInfo>::const_iterator p = registered_initial_conditions_models->begin();
           p != registered_initial_conditions_models->end(); ++p)
        std_cxx1x::get<1>(*p)(prm);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    template class Interface<deal_II_dimension>;

    template
    void
    register_initial_conditions_model<deal_II_dimension> (const std::string &,
                                                          void ( *) (ParameterHandler &),
                                                          Interface<deal_II_dimension> * ( *) ());

    template
    Interface<deal_II_dimension> *
    create_initial_conditions<deal_II_dimension> (ParameterHandler &prm,
                                                  const GeometryModel::Interface<deal_II_dimension> &geometry_model,
                                                  const BoundaryTemperature::Interface<deal_II_dimension> &boundary_temperature,
                                                  const AdiabaticConditions<deal_II_dimension>      &adiabatic_conditions);
  }
}
