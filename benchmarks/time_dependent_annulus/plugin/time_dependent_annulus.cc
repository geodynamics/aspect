/*
  Copyright (C) 2018 - 2019 by the authors of the ASPECT code.

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

#include "time_dependent_annulus.h"

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(TimeDependentAnnulus,
                                   "time dependent annulus",
                                   "This is the material model setup "
                                   "for the time dependent annular flow benchmark from Section 5 of "
                                   "``Gassmoeller, Lokavarapu, Bangerth, Puckett (2019): "
                                   "Evaluating the Accuracy of Hybrid Finite Element/Particle-In-Cell "
                                   "Methods for Modeling Incompressible Stokes Flow. Geophys. J. Int. "
                                   "submitted.'' "
                                   "The model "
                                   "can either use the analytic density, or whatever is stored in the "
                                   "first compositional field as the density that enters the Stokes equations. "
                                   "This allows testing different advection schemes for the density.")
  }

  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(TimeDependentAnnulus,
                                  "time dependent annulus",
                                  "This is the error postprocessor "
                                  "for the time dependent annular flow benchmark from Section 5 of "
                                  "``Gassmoeller, Lokavarapu, Bangerth, Puckett (2019): "
                                  "Evaluating the Accuracy of Hybrid Finite Element/Particle-In-Cell "
                                  "Methods for Modeling Incompressible Stokes Flow. Geophys. J. Int. "
                                  "submitted.'' "
                                  "Specifically, it can compute the errors for "
                                  "the velocity, pressure, and density variables.")
  }
}
