/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
/*  $Id$  */


#include <aspect/introspection.h>
#include <aspect/global.h>

namespace aspect
{

  template <int dim>
  Introspection<dim>::Introspection(const unsigned int n_compositional_fields)
    :
    extractors (n_compositional_fields)
  {
    // set up a mapping between vector components to the blocks they
    // correspond to. each variable has its own block except
    // for the velocities which are all mapped into block 0
    components_to_blocks.resize (dim+2+n_compositional_fields,0);
    components_to_blocks[dim] = 1;
    components_to_blocks[dim+1] = 2;
    for (unsigned int i=dim+2; i<dim+2+n_compositional_fields; ++i)
      components_to_blocks[i] = i-dim+1;
  }


  template <int dim>
  Introspection<dim>::Extractors::Extractors (const unsigned int n_compositional_fields)
    :
    velocities (0),
    pressure (dim),
    temperature (dim+1)
  {
    for (unsigned int c=0; c<n_compositional_fields; ++c)
      compositional_fields.push_back (FEValuesExtractors::Scalar(dim+2+c));
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Introspection<dim>;\
   
  ASPECT_INSTANTIATE(INSTANTIATE)
}
