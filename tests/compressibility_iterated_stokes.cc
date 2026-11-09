// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

// use the same postprocessing facilities as for the 'compressibility'
// testcase
#include "compressibility.cc"

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class CompressibilityIteratedStokes : public MaterialModel::Simple<dim>
    {
      public:

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
          * Return true if the compressibility() function returns something that
          * is not zero.
          */
        bool
        is_compressible () const override;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    CompressibilityIteratedStokes<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          out.densities[i] = 10.0/11.0*std::exp(in.pressure[i]/100.0);
          out.compressibilities[i] = 0.01;
        }
    }

    template <int dim>
    bool
    CompressibilityIteratedStokes<dim>::
    is_compressible () const
    {
      return true;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CompressibilityIteratedStokes,
                                   "compressibility iterated stokes",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a non-zero compressibility.")
  }
}
