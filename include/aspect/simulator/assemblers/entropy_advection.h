// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_simulator_assemblers_entropy_advection_h
#define _aspect_simulator_assemblers_entropy_advection_h


#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Assemblers
  {
    /**
     * This class assembles the terms for the matrix and right-hand-side of the entropy
     * advection equation for the current cell. For details on the purposes of the
     * individual functions of this class see Assemblers::Interface<dim>(). For
     * details on the entropy method see \cite dannberg:etal:2022 (https://doi.org/10.1093/gji/ggac293).
     */
    template <int dim>
    class EntropyAdvectionSystem : public Assemblers::Interface<dim>, public Assemblers::AdvectionStabilizationInterface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const override;

        std::vector<double>
        compute_residual(internal::Assembly::Scratch::ScratchBase<dim>  &scratch) const override;

        std::vector<double>
        advection_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const override;

        std::vector<double>
        diffusion_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const override;
    };

    /**
     * This function replaces the existing temperature assembler and the assembler
     * for the entropy composition equation with an assembler specifically designed
     * for the entropy method.
     */
    template <int dim>
    void set_assemblers_entropy_advection(const SimulatorAccess<dim> &simulator_access,
                                          Assemblers::Manager<dim> &assemblers);
  }
}


#endif
