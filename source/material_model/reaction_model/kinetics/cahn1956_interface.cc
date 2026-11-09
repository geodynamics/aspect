// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/material_model/reaction_model/kinetics/cahn1956_interface.h>
#include <aspect/global.h>
#include <deal.II/base/patterns.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      std::unique_ptr<Cahn1956Interface<dim>> create_reaction_model(const std::string &model_name)
      {
        return std::unique_ptr<Cahn1956Interface<dim>>(ReactionModelPluginList<dim>::create_plugin(model_name, "Reaction kinetics model"));
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
#define INSTANTIATE(dim) \
  template class Cahn1956Interface<dim>; \
  template std::unique_ptr<Cahn1956Interface<dim>> create_reaction_model<dim>(const std::string &);
      ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
    }
  }
}
