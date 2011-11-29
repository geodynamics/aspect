//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/geometry_model/interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

namespace aspect
{
  namespace GeometryModel
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


// -------------------------------- Deal with registering geometry models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      typedef
      std_cxx1x::tuple<std::string,
                std::string,
                void ( *) (ParameterHandler &),
                Interface<deal_II_dimension> * ( *) ()>
                GeometryModelInfo;

      // A pointer to a list of all postprocessors. the three elements of the tuple
      // correspond to the arguments given to the register_geometry_model
      // function.
      //
      // The object is a pointer rather for the same reason as discussed in
      // postprocess_base.cc for the corresponding variable there
      std::list<GeometryModelInfo> *registered_geometry_models = 0;
    }



    template <int dim>
    void
    register_geometry_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> * (*factory_function) ())
    {
      // see if this is the first time we get into this
      // function and if so initialize the variable above
      if (registered_geometry_models == 0)
        registered_geometry_models = new std::list<GeometryModelInfo>();

      // now add one record to the list
      registered_geometry_models->push_back (GeometryModelInfo(name,
                                                               description,
                                                               declare_parameters_function,
                                                               factory_function));
    }


    template <int dim>
    Interface<dim> *
    create_geometry_model (ParameterHandler &prm)
    {
      Assert (registered_geometry_models != 0, ExcInternalError());

      std::string model_name;
      prm.enter_subsection ("Geometry model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      for (std::list<GeometryModelInfo>::const_iterator p = registered_geometry_models->begin();
           p != registered_geometry_models->end(); ++p)
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
      Assert (registered_geometry_models != 0, ExcInternalError());

      // first collect a list of all registered models
      std::string model_names;
      for (std::list<GeometryModelInfo>::const_iterator p = registered_geometry_models->begin();
           p != registered_geometry_models->end(); ++p)
        {
          if (model_names.size() > 0)
            model_names += "|";
          model_names += std_cxx1x::get<0>(*p);
        }

      // then declare the actual entry in the parameter file
      prm.enter_subsection ("Geometry model");
      {
        prm.declare_entry ("Model name", "",
                           Patterns::Selection (model_names),
                           "Select one of the available geometry models.");
      }
      prm.leave_subsection ();

      for (std::list<GeometryModelInfo>::const_iterator p = registered_geometry_models->begin();
           p != registered_geometry_models->end(); ++p)
        std_cxx1x::get<2>(*p)(prm);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    template class Interface<deal_II_dimension>;

    template
    void
    register_geometry_model<deal_II_dimension> (const std::string &,
                                                const std::string &,
                                                void ( *) (ParameterHandler &),
                                                Interface<deal_II_dimension> * ( *) ());

    template
    Interface<deal_II_dimension> *
    create_geometry_model<deal_II_dimension> (ParameterHandler &prm);
  }
}
