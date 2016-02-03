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


#ifndef __aspect__melt_h
#define __aspect__melt_h

#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/assembly.h>

namespace aspect
{
  using namespace dealii;


  namespace Assemblers
  {
    /**
     * A class for the definition of functions that implement the
     * linear system terms for the *melt* migration compressible or
     * incompressible equations.
     */
    template <int dim>
    class MeltEquations : public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>,
      public SimulatorAccess<dim>
    {
      public:

        /**
         * Compute the integrals for the preconditioner for the Stokes system in
         * the case of melt migration on a single cell.
         */
        void
        local_assemble_stokes_preconditioner_melt (const double                                             pressure_scaling,
                                                   internal::Assembly::Scratch::StokesPreconditioner<dim>  &scratch,
                                                   internal::Assembly::CopyData::StokesPreconditioner<dim> &data) const;

        /**
         * Compute the integrals for the Stokes matrix and right hand side in
         * the case of melt migration on a single cell.
         */
        void
        local_assemble_stokes_system_melt (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           const double                                     pressure_scaling,
                                           const bool                                       rebuild_stokes_matrix,
                                           internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                           internal::Assembly::CopyData::StokesSystem<dim> &data) const;

        /**
         * Compute the boundary integrals for the Stokes right hand side in
         * the case of melt migration on a single cell. These boundary terms
         * are used to describe Neumann boundary conditions for the fluid pressure.
         */
        void
        local_assemble_stokes_system_melt_boundary (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                    const unsigned int                                    face_no,
                                                    const double                                          pressure_scaling,
                                                    internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                                    internal::Assembly::CopyData::StokesSystem<dim>      &data) const;

      private:
        /**
         * Returns the right hand side of the fluid pressure equation in
         * the case of melt migration for a single quadrature point.
         */
        double
        compute_fluid_pressure_RHS(const internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                   MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                   MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                                   const unsigned int q_point) const;
    };

    /**
     * A namespace for the definition of functions that implement various
     * other terms that need to occasionally or always be assembled.
     */
    namespace OtherTerms
    {
      /**
       * Integrate the local fluid pressure shape functions on a single cell
       * for models with melt migration, so that they can later be used to do
       * the pressure right-hand side compatibility modification.
       */
      template <int dim>
      void
      pressure_rhs_compatibility_modification_melt (const SimulatorAccess<dim>                      &simulator_access,
                                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                    internal::Assembly::CopyData::StokesSystem<dim> &data);
    }
  }
}

#endif
