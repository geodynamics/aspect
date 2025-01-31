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

#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/material_model/simpler.h>

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    class CellMaterial :
      public aspect::MaterialModel::Simpler<dim>, aspect::SimulatorAccess<dim>
    {
      public:
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          Simpler<dim>::evaluate(in,out);

          if (in.current_cell.state() == IteratorState::valid)
            {
              std::cout << "Level: " << in.current_cell->level() << " Index: " << in.current_cell->index() << std::endl;
            }
        }
    };
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CellMaterial,
                                   "cell",
                                   "")
  }
}
