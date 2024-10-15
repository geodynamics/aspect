/*
  Copyright (C) 2017 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_simulator_assemblers_advection_h
#define _aspect_simulator_assemblers_advection_h


#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Assemblers
  {
    /**
     * This class assembles the terms for the matrix and right-hand-side of the advection
     * equation for the current cell.
     */
    template <int dim>
    class AdvectionSystem : public Assemblers::Interface<dim>, public Assemblers::AdvectionStabilizationInterface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        std::vector<double>
        compute_residual(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base) const override;
    };

    template <int dim>
    class DarcySystem : public Assemblers::Interface<dim>, public Assemblers::AdvectionStabilizationInterface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        std::vector<double>
        compute_residual(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base) const override;

        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side equation for the
     * current cell in case we only want to solve the diffusion equation.
     */
    template <int dim>
    class DiffusionSystem : public Assemblers::Interface<dim>, public Assemblers::AdvectionStabilizationInterface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        std::vector<double>
        compute_residual(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base) const override;

        virtual
        std::vector<double>
        advection_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const override;

        virtual
        std::vector<double>
        diffusion_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const override;
    };

    /**
     * This class assembles the face terms for the right-hand-side of the
     * advection equation for a face at the boundary of the domain where
     * Neumann boundary conditions are used (which allow to prescribe a heat flux).
     */
    template <int dim>
    class AdvectionSystemBoundaryHeatFlux : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the face terms for the matrix and right-hand-side of
     * the discontinuous advection equation for a face at the boundary of the domain.
     */
    template <int dim>
    class AdvectionSystemBoundaryFace : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the face terms for the matrix and right-hand-side of
     * the discontinuous advection equation for a face in the interior of the domain.
     */
    template <int dim>
    class AdvectionSystemInteriorFace : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };
  }
}


#endif
