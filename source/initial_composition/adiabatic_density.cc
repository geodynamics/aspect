/*
  Copyright (C) 2017 - 2023 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/initial_composition/adiabatic_density.h>

#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    double
    AdiabaticDensity<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      if (n_comp == this->introspection().find_composition_type(CompositionalFieldDescription::density))
        return this->get_adiabatic_conditions().density(position);

      return 0.0;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(AdiabaticDensity,
                                              "adiabatic density",
                                              "Specify the initial composition as the adiabatic reference density at "
                                              "each position. Note that only the field of the type 'density' "
                                              "will be filled. For all other fields this plugin returns 0.0.")
  }
}
