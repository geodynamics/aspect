//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/model.h>

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

}

// explicit instantiations
namespace aspect
{

  template class MaterialModel<deal_II_dimension>;

}
