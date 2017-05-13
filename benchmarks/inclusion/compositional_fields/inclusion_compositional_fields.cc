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
