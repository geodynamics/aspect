// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
     * This class assembles the face terms for the left-hand and right-hand side of the
     * advection equation for a face at the boundary of the domain where
     * Robin boundary conditions are used (which allow to prescribe a heat flux and a temperature).
     */
    template <int dim>
    class AdvectionSystemRobinBoundary : public Assemblers::Interface<dim>,
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
