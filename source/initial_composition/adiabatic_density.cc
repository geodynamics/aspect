// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
      if (n_comp == this->introspection().find_composition_type(CompositionalFieldDescription::density) &&
          this->get_adiabatic_conditions().is_initialized())
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
