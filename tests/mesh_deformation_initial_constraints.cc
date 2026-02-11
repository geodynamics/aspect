/*
  Copyright (C) 2020 - 2024 by the authors of the ASPECT code.

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
#include <deal.II/dofs/dof_tools.h>

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

        virtual
        void
        compute_initial_deformation_as_constraints(const Mapping<dim> &/*mapping*/,
                                                   const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                   const types::boundary_id boundary_indicator,
                                                   AffineConstraints<double> &constraints) const override
        {

          const IndexSet constrained_dofs =
            DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                            ComponentMask(dim, true),
          {boundary_indicator});

          for (const types::global_dof_index index : constrained_dofs)
            {
              if (constraints.can_store_line(index))
                if (constraints.is_constrained(index)==false)
                  {
                    // set some nonsensical values:
                    const double value = static_cast<double>(index) * 10.0;
#if DEAL_II_VERSION_GTE(9,6,0)
                    constraints.add_constraint(index,
                                               {},
                                               value);
#else
                    constraints.add_line(index);
                    constraints.set_inhomogeneity(index, value);
#endif
                  }
            }

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
