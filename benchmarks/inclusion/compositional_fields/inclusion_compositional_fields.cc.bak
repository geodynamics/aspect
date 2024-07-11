/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
#include "inclusion_compositional_fields.h"



// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(InclusionCompositionalMaterial,
                                   "InclusionCompositionalMaterial",
                                   "A material model that corresponds to the 'Inclusion' benchmark "
                                   "defined in Duretz et al., G-Cubed, 2011. using "
                                   "compositional fields.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(InclusionBoundary,
                                            "InclusionBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "``inclusion'' benchmark. See the manual and the Kronbichler, Heister "
                                            "and Bangerth paper on ASPECT for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(InclusionPostprocessor,
                                  "InclusionPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Duretz et al., G-Cubed, 2011, paper with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")
  }
}
