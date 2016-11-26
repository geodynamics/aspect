/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__assemblers_h
#define __aspect__assemblers_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  template <int dim>
  class AdvectionAssembler : public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>,
    public SimulatorAccess<dim>
  {
    public:

      void
      local_assemble_advection_system (const typename Simulator<dim>::AdvectionField &advection_field,
                                       const double artificial_viscosity,
                                       internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                       internal::Assembly::CopyData::AdvectionSystem<dim> &data) const;

      std::vector<double>
      compute_advection_system_residual(const typename Simulator<dim>::AdvectionField     &advection_field,
                                        internal::Assembly::Scratch::AdvectionSystem<dim> &scratch) const;

      void
      local_assemble_discontinuous_advection_boundary_face_terms(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                 const unsigned int face_no,
                                                                 const typename Simulator<dim>::AdvectionField &advection_field,
                                                                 internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                                                 internal::Assembly::CopyData::AdvectionSystem<dim> &data) const;

      void
      local_assemble_discontinuous_advection_interior_face_terms(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                 const unsigned int face_no,
                                                                 const typename Simulator<dim>::AdvectionField &advection_field,
                                                                 internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                                                 internal::Assembly::CopyData::AdvectionSystem<dim> &data) const;
  };


  template <int dim>
  class StokesAssembler : public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>,
    public SimulatorAccess<dim>
  {
    public:
      void
      local_assemble_stokes_preconditioner (const double                                             pressure_scaling,
                                            internal::Assembly::Scratch::StokesPreconditioner<dim>  &scratch,
                                            internal::Assembly::CopyData::StokesPreconditioner<dim> &data) const;
      void
      local_assemble_stokes_incompressible (const double                                     pressure_scaling,
                                            const bool                                       rebuild_stokes_matrix,
                                            internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                            internal::Assembly::CopyData::StokesSystem<dim> &data) const;
      void
      local_assemble_stokes_compressible_diffusion (const double                                     pressure_scaling,
                                                    const bool                                       rebuild_stokes_matrix,
                                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                    internal::Assembly::CopyData::StokesSystem<dim> &data) const;

      void
      local_assemble_stokes_mass_density_gradient (const double                                     pressure_scaling,
                                                   const bool                                       rebuild_stokes_matrix,
                                                   internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                   internal::Assembly::CopyData::StokesSystem<dim> &data,
                                                   const Parameters<dim> &parameters) const;

      void
      local_assemble_stokes_mass_density_gradient_implicit (const double                                     pressure_scaling,
                                                            const bool                                       rebuild_stokes_matrix,
                                                            internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                            internal::Assembly::CopyData::StokesSystem<dim> &data,
                                                            const Parameters<dim> &parameters) const;

      void
      local_assemble_stokes_mass_density_explicit (const double                                     pressure_scaling,
                                                   const bool                                       rebuild_stokes_matrix,
                                                   internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                   internal::Assembly::CopyData::StokesSystem<dim> &data,
                                                   const Parameters<dim> &parameters) const;

  };
}

#endif
