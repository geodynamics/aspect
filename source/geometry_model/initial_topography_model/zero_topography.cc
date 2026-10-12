// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2017 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/geometry_model/initial_topography_model/zero_topography.h>


namespace aspect
{
  namespace InitialTopographyModel
  {
    template <int dim>
    double
    ZeroTopography<dim>::
    value (const Point<dim-1> &/*p*/) const
    {
      // return a zero value regardless of position
      return 0.0;
    }

    template <int dim>
    double
    ZeroTopography<dim>::
    max_topography () const
    {
      return 0;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(ZeroTopography,
                                             "zero topography",
                                             "Implementation of a model in which the initial topography "
                                             "is zero. ")
  }
}
