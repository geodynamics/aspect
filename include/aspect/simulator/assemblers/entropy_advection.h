/*
  Copyright (C) 2016 - 2025 by the authors of the ASPECT code.

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
