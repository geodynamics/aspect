/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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
