//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/model_base.h>
#include <aspect/model_simple.h>
#include <aspect/model_table.h>
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


    template <int dim>
    Interface<dim> *
    create (const std::string &name)
    {
      if (name == "simple")
        return new Simple<dim>();
      else if (name == "table")
        return new Table<dim>();
      else
        throw dealii::ExcNotImplemented();
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
    Interface<deal_II_dimension> *
    create<deal_II_dimension> (const std::string &);
  }
}
