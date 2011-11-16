//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/material_model_base.h>
#include <deal.II/base/exceptions.h>

namespace aspect
{
  namespace MaterialModel
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


// -------------------------------- Deal with registering material models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      typedef
      std_cxx1x::tuple<std::string,
                void ( *) (ParameterHandler &),
                Interface<deal_II_dimension> * ( *) ()>
                MaterialModelInfo;

      // A pointer to a list of all postprocessors. the three elements of the tuple
      // correspond to the arguments given to the register_material_model
      // function.
      //
      // The object is a pointer rather for the same reason as discussed in
      // postprocess_base.cc for the corresponding variable there
      std::list<MaterialModelInfo> *registered_material_models = 0;
    }



    template <int dim>
    void
    register_material_model (const std::string &name,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> * (*factory_function) ())
    {
      // see if this is the first time we get into this
      // function and if so initialize the variable above
      if (registered_material_models == 0)
        registered_material_models = new std::list<MaterialModelInfo>();

      // now add one record to the list
      registered_material_models->push_back (MaterialModelInfo(name,
                                                               declare_parameters_function,
                                                               factory_function));
    }


    template <int dim>
    Interface<dim> *
    create_material_model (const std::string &name)
    {
      Assert (registered_material_models != 0, ExcInternalError());

      for (std::list<MaterialModelInfo>::const_iterator p = registered_material_models->begin();
           p != registered_material_models->end(); ++p)
        if (std_cxx1x::get<0>(*p) == name)
          return std_cxx1x::get<2>(*p)();

      AssertThrow (false, ExcNotImplemented());
      return 0;
    }



    void
    declare_parameters (ParameterHandler &prm)
    {
      Assert (registered_material_models != 0, ExcInternalError());

      for (std::list<MaterialModelInfo>::const_iterator p = registered_material_models->begin();
           p != registered_material_models->end(); ++p)
        std_cxx1x::get<1>(*p)(prm);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    template class Interface<deal_II_dimension>;

    template
    void
    register_material_model<deal_II_dimension> (const std::string &,
                                                void ( *) (ParameterHandler &),
                                                Interface<deal_II_dimension> * ( *) ());

    template
    Interface<deal_II_dimension> *
    create_material_model<deal_II_dimension> (const std::string &);
  }
}
