#include "solkz.h"



// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolKzMaterial,
                                   "SolKzMaterial",
                                   "A material model that corresponds to the 'SolKz' benchmark "
                                   "defined in Duretz et al., G-Cubed, 2011.")


    ASPECT_REGISTER_POSTPROCESSOR(SolKzPostprocessor,
                                  "SolKzPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Duretz et al., G-Cubed, 2011, paper with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")
  }
}
