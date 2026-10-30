// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2020 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
