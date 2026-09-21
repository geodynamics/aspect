// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/mesh_deformation/interface.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    class TestMeshDeformation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void save (std::map<std::string, std::string> &status_strings) const override
        {
          this->get_pcout() << "checkpoint_id: " << this->get_checkpoint_id() << std::endl;
        }
    };
  }
}


namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(TestMeshDeformation,
                                           "Test mesh deformation",
                                           "Test mesh deformation plugin for checking that the checkpoint id is equal to "
                                           "the current checkpoint id in the simulator when the save() function is called.");
  }
}
