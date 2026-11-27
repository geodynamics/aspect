/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#include <aspect/material_model/reaction_model/gibbs_free_energy/interface.h>
#include <aspect/global.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace GibbsFreeEnergyModel
      {
        template <int dim>
        std::unique_ptr<Interface<dim>> create_gibbs_free_energy_model(const std::string &model_name)
        {
          return std::unique_ptr<Interface<dim>>(PluginList<dim>::create_plugin(model_name, "Gibbs free energy model"));
        }
      }
    }
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace GibbsFreeEnergyModel
      {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template std::unique_ptr<Interface<dim>> create_gibbs_free_energy_model<dim>(const std::string &);
        ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
      }
    }
  }
}
