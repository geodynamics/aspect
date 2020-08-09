/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {
    template<int dim>
    class PrescribedDeformation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        PrescribedDeformation()
        {}

        Tensor<1,dim>
        compute_initial_deformation_on_boundary(const types::boundary_id /*boundary_indicator*/,
                                                const Point<dim> &position) const
        {
          const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(position);
          Tensor<1,dim> topography_direction;
          if (gravity.norm() > 0.0)
            topography_direction = -gravity / gravity.norm();

          const double topography_amplitude = 20000. * std::cos(2.*numbers::PI*position[0]/(660000.));
          return topography_amplitude * topography_direction;
        }
    };
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(PrescribedDeformation,
                                           "prescribed deformation",
                                           "A test plugin for initial mesh deformation.")
  }
}
