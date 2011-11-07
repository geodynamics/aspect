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
}

// explicit instantiations
namespace aspect
{

  template class MaterialModel<deal_II_dimension>;

}
