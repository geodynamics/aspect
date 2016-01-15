/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/melt_interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    MeltInterface<dim>::~MeltInterface ()
    {}



    template <int dim>
    MeltMaterialModelOutputs<dim>::MeltMaterialModelOutputs(const unsigned int n_points,
                                                            const unsigned int n_comp)
      : MaterialModelOutputs<dim>(n_points, n_comp)
    {
      compaction_viscosities.resize(n_points);
      fluid_viscosities.resize(n_points);
      permeabilities.resize(n_points);
      fluid_densities.resize(n_points);
      fluid_compressibilities.resize(n_points);
    }

  }
}

// explicit instantiations
namespace aspect
{

  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class MeltInterface<dim>;\
  template class MeltMaterialModelOutputs<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
