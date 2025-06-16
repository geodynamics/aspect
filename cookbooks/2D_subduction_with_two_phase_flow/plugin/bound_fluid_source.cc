/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

#include <aspect/material_model/reactive_fluid_transport.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class BoundFluidSource: public MaterialModel::ReactiveFluidTransport<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const override;
    };


    template <int dim>
    void
    BoundFluidSource<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      ReactiveFluidTransport<dim>::evaluate(in, out);
      const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid_idx");
      for (unsigned int q=0; q < in.n_evaluation_points(); ++q)
        if (out.reaction_terms[q][bound_fluid_idx] <= 0.0)
          out.reaction_terms[q][bound_fluid_idx] = 0.0;
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(BoundFluidSource,
                                   "reactive fluid transport bound fluid source",
                                   "A simple material model that is like the "
                                   "'reactive fluid transport' model, but prevents "
                                   "the bound fluid content from decreasing.")
  }
}
