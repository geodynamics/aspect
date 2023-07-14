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

#include <aspect/material_model/interface.h>
#include <aspect/global.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  template <int dim>
  class TestMaterial:
    public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            out.viscosities[i] = 1e19;
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1e-4;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 3000;
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              {
                if (in.requests_property(MaterialModel::MaterialProperties::reaction_terms))
                  out.reaction_terms[i][c] = trace(in.strain_rate[i]) * this->get_timestep();
                else
                  out.reaction_terms[i][c] = 0.0;
              }

          }
      }
  };
}

// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(TestMaterial,
                                 "test material",
                                 "")

}
