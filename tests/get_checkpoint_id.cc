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
