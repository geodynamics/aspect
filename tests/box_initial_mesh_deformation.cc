// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/mesh_deformation/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
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
