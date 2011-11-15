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
  template <int dim>
  MaterialModel<dim>::~MaterialModel ()
  {}


  template <int dim>

  void
  MaterialModel<dim>::
  declare_parameters (dealii::ParameterHandler &prm)
  {}


  template <int dim>
  void
  MaterialModel<dim>::parse_parameters (dealii::ParameterHandler &prm)
  {}

  template <int dim>
  MaterialModel<dim> *
  MaterialModel<dim>::create( const std::string name)
  {
    if (name == "simple")
      return new MaterialModel_Simple<dim>();
    else if (name == "table")
      return new MaterialModel_Table<dim>();
    else
      throw dealii::ExcNotImplemented();
  }

}

// explicit instantiations
namespace aspect
{

  template class MaterialModel<deal_II_dimension>;

}
