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


#include <aspect/initial_conditions/box.h>
#include <aspect/geometry_model/box.h>


namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    double
    PerturbedBox<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // spherical shell. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (this->geometry_model);
      Assert (geometry != 0,
              ExcMessage ("This initial condition can only be used if the geometry "
                          "is a box."));

      double perturbation = 1;
      for (unsigned int d=0; d<dim; ++d)
        perturbation *= std::sin(numbers::PI*position[d]/geometry->get_extents()[d]);
      return 1 + perturbation/10;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    template class PerturbedBox<deal_II_dimension>;
    ASPECT_REGISTER_INITIAL_CONDITIONS(PerturbedBox,
                                       "perturbed box",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed slightly from an otherwise constant value "
                                       "equal to one. The perturbation is chosen in such a way "
                                       "that the initial temperature is constant to one along "
                                       "the entire boundary.");
  }
}
